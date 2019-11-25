# Scotch

Scotch is a machine-learning caller for indels.

## Installation 

Clone this repository. And clone the repository with pre-calculated features that describe the reference genome for the genome build you're using. 

```
$ git clone https://github.com/AshleyLab/scotch.git
$ git clone https://github.com/AshleyLab/scotch-data-grch37 # for GRCh37
$ git clone https://github.com/AshleyLab/scotch-data-grch38 # for GRCh38
```

## Dependencies

Scotch requires the following command-line utilities:
```
samtools
```

Scotch requires the following Python packages:

```
numpy
pysam
```
Scotch requires the following R packages:
```
data.table
reshape
randomForest
```

## Run 

### Overview

```
# extract downloaded region features for regions of interest
# ~/scotch-data/ is clone of https://github.com/AshleyLab/scotch-data-grch37 or https://github.com/AshleyLab/scotch-data-grch38
python ~/scotch/scotch.py prepare-region-features --beds_dir=~/beds/ --all_rfs_dir=~/scotch-data/ --output_trim_rfs_dir=~/trim_rfs/

# remove duplicate reads from input BAM
python ~/scotch/scotch.py rmdup-bam --project_dir=~/ABC123/ --bam=~/input.bam

# can process BAM in parallel by chromosome
for chr in {1..22} X Y
do
	# create bam with soft clipping reverted
	python ~/scotch/scotch.py unclip-bam --project_dir=~/ABC123/ --chrom=$chr --beds_dir=~/beds/ --fasta_ref=~/GRCh37.fa
	
	# can calculate features in parallel	
	for feature in {get-features-depth,get-features-nReads,get-features-read}
	do
		python ~/scotch/scotch.py $feature --project_dir=~/ABC123/ --chrom=$chr --beds_dir=~/beds/ --fasta_ref=~/GRCh37.fa
	done
	
	# compile features
	python ~/scotch/scotch.py compile-features --project_dir=~/ABC123/ --chrom=$chr --beds_dir=~/beds/ --trim_rfs_dir=~/trim_rfs/
	
	# make predictions
	python ~/scotch/scotch.py predict --project_dir=~/ABC123/ --chrom=$chr --fasta_ref=~/GRCh37.fa
done

```

### Input
Scotch accepts a Binary Alignment Mapping (BAM) file containing whole-genome next generation sequencing data. Scotch also accepts a FASTA file providing the corresponding reference genome. Scotch divides the input by chromosome for parallel processing. 

### Common arguments

* `--project_dir`: for a run of Scotch, where intermediate and output files should be stored
* `--beds_dir`: a directory of BED files  specifying the genomic regions Scotch should analyze
  * One file for each chrom `{{1..22},X,Y}.bed` (e.g., `4.bed`)
  * Where each line of each file is in the format `chrom\tstart\tstop` to indicate the region `chrom:start-stop` (e.g., `4:12822-1423146`)
    * Note: the expected format of `chrom` is `4` without the "chr", not `chr4`

## Pipeline stages

### `prepare-region-features`
##### Required args
* `--beds_dir`
* `--all_rfs_dir`: path to directory with all region features, probably location where `https://github.com/AshleyLab/scotch-data-grch37` or `https://github.com/AshleyLab/scotch-data-grch38` was cloned
* `--output_trim_rfs_dir`: path to directory where should output trimmed region features; can be an empty directory

Scotch's model relies on eight *region features* that describe the reference genome. For convenience, we've computed them across the entire GRCh37 reference genome and made them available in two other repositories. Run this command with the bed files that describe your regions of interest to obtain the corresponding region features. 

### `rmdup-bam`
##### Required args
* `--project_dir`
* `--bam`: path to indexed, BWA-aligned, whole-genome, next-generation input sequencing data

Scotch removes duplicate reads from the input bam with `samtools rmdup`.

### `unclip-bam`
##### Required args
* `--project_dir`
* `--chrom`
* `--beds_dir`
* `--fasta_ref`
##### Optional args
* `--gatk_jar` (default `.../aux/gatk-3.8.jar`): Scotch requires GATK 3.8 (it uses commands that have not yet been ported over to GATK 4.0) and is shipped with a version of this, but if you'd like to use your own, the path to the `jar` can be specified here
* `--gatk_mem` (default: `5`): the amount of memory, in GB, to run GATK with

Many of Scotch's features relate to soft clipping. This state of the pipeline creates, for a single chromosome, a BAM from the output of `rmdup-bam`, but with all soft clipped bases reverted. (This stage's output is used in `get-features-nReads` to calculate the coverage including soft clipping, whereas `get-features-depth` calculates coverage excludings soft clipping). 

### `get-features-depth`
##### Required args
* `--project_dir`
* `--chrom`
* `--beds_dir`
* `--fasta_ref` 
##### Optional args
* `--gatk_jar` (default `.../aux/gatk-3.8.jar`)
* `--gatk_mem` (default: `5`)

Calculates coverage, within the specified chromosome, excluding soft clipping. In the directory `{project_dir}/features/{chrom}/`, produces `depth.feat.gz`, `depth.feat.log`, and `depth.feat.stats`. 

### `get-features-nReads`
##### Required args
* `--project_dir`
* `--chrom`
* `--beds_dir`
* `--fasta_ref` 
##### Optional args
* `--gatk_jar` (default `.../aux/gatk-3.8.jar`)
* `--gatk_mem` (default: `5`)

Calculates coverage, within the specified chromosome, including soft clipping. In the directory `{project_dir}/features/{chrom}/`, produces `nReads.feat.gz`, `nReads.feat.log`, and `nReads.feat.stats`. 

### `get-features-read`
##### Required args
* `--project_dir`
* `--chrom`
* `--beds_dir`
* `--fasta_ref` 

Calculates several features within the specified chromosome, including (per-position) mean mapping quality, mean base quality, operations described in reads' CIGAR strings, and several features describing soft clipping. In the directory `{project_dir}/features/{chrom}/`, produces `read.feats.gz`.


### `compile-features`
##### Required args
* `--project_dir`
* `--chrom`
* `--beds_dir`
* `--trim_rfs_dir`: path to directory with trimmed region features, produced in `prepare-region-features` as `--output_trim_rfs_dir`

Combines the results of `get-features-depth`, `get-features-nReads`, and `get-features-read`. Normalizes some features for inter-sample comparability, and computes new features based on combinations of distinct features.  In the directory `{project_dir}/features/{chrom}/`, produces `matrix.txt.gz`.


### `predict`
* `--project_dir`
* `--chrom`
* `--fasta_ref` 

From the feature matrix, makes predictions against Scotch's random forest model. In the directory `{project_dir}/results/`, produces `results.{chrom}.vcf` and other output files described further below. 


### Output

Scotch produces several output files. `scotch.{chrom}.vcf` includes all the results in VCF format. Alternate alleles may be represented as`<DEL_L>`,`<DEL_R>` or `<INS>` representing a deletion start, deletion end or insertion breakpoint, respectively. 

#### Encoding

Scotch also produces several VCF files where indel breakpoints are _encoded_ as regular variants. The motivation is that some tools (including Scotch and Pindel) do not report the nucleotide sequence of alternate alleles for all variants. (They may, for example, report just `<INS>` instead.) As a result, output VCFs that include their calls may not be recognized as valid VCFs. 

`encode.py` translates each indel breakpoint into an SNV at the same locus with an arbitary alternate allele, producing strictly valid VCF output. `{stub}.encode_del_L.vcf` includes deletion start breakpoints represented this way, `{stub}.encode_del_R.vcf` includes deletion end breakpoints, `{stub}.encode_ins.vcf` includes insertion breakpoints, and `{stub}.encode_all.vcf` includes all breakpoints. 

Since this process preserves breakpoint position, these files can be input to benchmarking tools like GA4GH Benchmarking on precisionFDA, which execute a distance-based comparison to evaluate tools' performance. Truth VCFs and the VCFs output by other callers to be benchmarked should also be encoded by `encode.py`. The script is called as

```
python encode.py input.vcf output_stub reference.fa
```

## Features

Scotch’s model evaluates each position with respect to 40 features. These include “primary metrics,” quantities which are extracted directly from sequencing data; “delta features” which track the differences in primary features between neighboring positions; and “genomic features,” which describe the content of the reference genome at a given locus. 

#### Primary features
These 12 features are calculated directly from the sequencing data. Three describe coverage—including the number of reads, reads with no soft-clipping, and reads with a base quality of 13 or higher. Each of these are normalized across the sample for comparability with samples from various sequencing runs. Two more features describe the quality of the sequencing—the mean base quality and the mean mapping quality across all reads. Four more are calculated from the CIGAR string that details each read’s alignment to the reference—recording the proportion of bases at that position across all reads that are marked as inserted, deleted, soft-clipped, and that at are at the boundary of soft-clipping (i.e., the base is soft-clipped but at least one neighboring base is not). Two more features describe the soft-clipping of the reads, if present: one gives the mean base quality of soft-clipped bases, another gives the *consistency score* of the soft-clipping. 

A position’s consistency score is a metric we derived that gives the ratio of the number of reads supporting the most common soft-clipped base (i.e., A, T, C, or G), to the number of all soft-clipped reads. Soft-clipping provides important signal of an indel to our model; this score helps a model distinguish indel-related soft-clipping (where all soft-clipped reads should support the same nucleotide) from that caused by low sequencing quality (where different nucleotides will be present). 

#### Delta features
20 additional features give the change in each of the primary features listed above— except the soft-clipping consistency score—from a given locus to both of its neighbors.  

#### Genomic features
Eight features, lastly, are derived from the reference genome, providing Scotch with insight into regions where sequencing errors are more common. Four of these features are binary: they indicate whether a genomic position is located in high-confidence regions, “superdup” regions, repetitive regions, and low-complexity regions. The remaining four describe GC-content (in windows of 50 and 1000 bp), mappability, and uniqueness. For GRCh37, they're available at `https://github.com/AshleyLab/scotch-data-grch37`. For GRCh38, they're availalbe at `https://github.com/AshleyLab/scotch-data-grch38`.
