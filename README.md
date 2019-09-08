# Scotch

## Installation 

Clone this repository. And clone the [repository](https://github.com/AshleyLab/scotch-data) with pre-calculated features that describe the reference genome. 

```
$ git clone https://github.com/AshleyLab/scotch.git
$ git clone https://github.com/AshleyLab/scotch-data
```
## Run 

### Overview

```
# extract downloaded region features for GRCh37 regions of interest
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

*GRCH37*??

#### Output

### Common arguments

* `--project_dir` (exc. for prepare-region-features)

### Pipeline stages

#### `prepare-region-features`
##### Required args
* `--beds_dir`
* `--all_rfs_dir`: path to directory with all region features, probably location where `https://github.com/AshleyLab/scotch-data` was cloned
* `--output_trim_rfs_dir`: path to directory where should output trimmed region features; can be an empty directory

Scotch's model relies on eight *region features*. Unlike other features that describe the input sample, these features desribe the reference genome. For example, they include GC content, uniqueness, and mappability. Since region features are the same for all samples, we've computed them across the entire GRCh37 reference genome and made them available in another [repository](https://github.com/AshleyLab/scotch-data). But a user might not be interested in the entire genome. So given a directory of bed files that describe the regions of interest, Scotch subsets these portions of the region features and outputs them to another directory for use later, in `compile-features`. 

#### `rmdup-bam`
##### Required args
* `--project_dir`
* `--bam`: path to indexed, BWA-aligned, whole-genome, next-generation input sequencing data

Scotch removes duplicate reads from the input bam with `samtools rmdup`. We found the best results with calculating our features from such a bam with duplicates removed. 

#### `unclip-bam`
##### Required args
* `--project_dir`
* `--chrom`
* `--beds_dir`
* `--fasta_ref`
##### Optional args
* `--gatk_jar` (default `.../aux/gatk-3.8.jar`): Scotch requires GATK 3.8 (it uses commands that have not yet been ported over to GATK 4.0) and is shipped with a version of this, but if you'd like to use your own, the path to the `jar` can be specified here
* `--gatk_mem` (default: `5`): the amount of memory, in GB, to run GATK with

Many of Scotch's features relate to soft clipping. This state of the pipeline creates, for a single chromosome, a BAM from the output of `rmdup-bam`, but with all soft clipped bases reverted. (This stage's output is used in `get-features-nReads` to calculate the coverage including soft clipping, whereas `get-features-depth` calculates coverage excludings soft clipping). 

#### `get-features-depth`
##### Required args
* `--project_dir`
* `--chrom`
* `--beds_dir`
* `--fasta_ref` 
##### Optional args
* `--gatk_jar` (default `.../aux/gatk-3.8.jar`)
* `--gatk_mem` (default: `5`)

Calculates coverage, within the specified chromosome, excluding soft clipping. In the directory `{project_dir}/features/{chrom}/`, produces `depth.feat.gz`, `depth.feat.log`, and `depth.feat.stats`. 

#### `get-features-nReads`
##### Required args
* `--project_dir`
* `--chrom`
* `--beds_dir`
* `--fasta_ref` 
##### Optional args
* `--gatk_jar` (default `.../aux/gatk-3.8.jar`)
* `--gatk_mem` (default: `5`)

Calculates coverage, within the specified chromosome, including soft clipping. In the directory `{project_dir}/features/{chrom}/`, produces `nReads.feat.gz`, `nReads.feat.log`, and `nReads.feat.stats`. 

#### `get-features-read`
##### Required args
* `--project_dir`
* `--chrom`
* `--beds_dir`
* `--fasta_ref` 

Calculates several features within the specified chromosome, including (per-position) mean mapping quality, mean base quality, operations described in reads' CIGAR strings, and several features describing soft clipping. In the directory `{project_dir}/features/{chrom}/`, produces `read.feats.gz`.


#### `compile-features`
##### Required args
* `--project_dir`
* `--chrom`
* `--beds_dir`
* `--trim_rfs_dir`: path to directory with trimmed region features, produced in `prepare-region-features` as `--output_trim_rfs_dir`

Combines the results of `get-features-depth`, `get-features-nReads`, and `get-features-read`. Normalizes some features for inter-sample comparability, and computes new features based on combinations of distinct features.  In the directory `{project_dir}/features/{chrom}/`, produces `matrix.txt.gz`.


#### `predict`
* `--project_dir`
* `--chrom`
* `--fasta_ref` 

From the feature matrix, makes predictions against Scotch's random forest model. 

## More

### Features
Scotch’s model evaluates each position with respect to 39 features. These include “primary metrics,” quantities which are extracted directly from sequencing data; “delta features” which track the differences in primary features between neighboring positions; and “genomic features,” which describe the content of the reference genome at a given locus. Information on feature importance is available in the Supplementary Note (Supplementary Fig. 1, Supplementary Table 25). 

#### Primary features
These 11 features are calculated directly from the sequencing data. Three describe coverage—including the number of reads, reads with no soft-clipping, and reads with a base quality of 13 or higher. Each of these are normalized across the sample for comparability with samples from various sequencing runs. Two more features describe the quality of the sequencing—the mean base quality and the mean mapping quality across all reads. Four more are calculated from the CIGAR string that details each read’s alignment to the reference—recording the proportion of bases at that position across all reads that are marked as inserted, deleted, soft-clipped, and that at are at the boundary of soft-clipping (i.e., the base is soft-clipped but at least one neighboring base is not). Two more features describe the soft-clipping of the reads, if present: one gives the mean base quality of soft-clipped bases, another gives the *consistency score* of the soft-clipping. 

A position’s consistency score is a metric we derived that gives the ratio of the number of reads supporting the most common soft-clipped base (i.e., A, T, C, or G), to the number of all soft-clipped reads. Soft-clipping provides important signal of an indel to our model; this score helps a model distinguish indel-related soft-clipping (where all soft-clipped reads should support the same nucleotide) from that caused by low sequencing quality (where different nucleotides will be present). 

#### Delta features
20 additional features give the change in each of the primary features listed above— except the soft-clipping consistency score—from a given locus to both of its neighbors.  

#### Genomic features
Eight features, lastly, are derived from the reference genome, providing Scotch with insight into regions where sequencing errors are more common. Four of these features are binary: they indicate whether a genomic position is located in high-confidence regions, “superdup” regions, repetitive regions, and low-complexity regions. The remaining four describe GC-content (in windows of 50 and 1000 bp), mappability, and uniqueness. 

### Prediction and Output
These features are combined in a human-readable TSV that can serve as the input to any number of machine-learning setups. We trained several random forest models to identify the signals of indels in this data. The primary output of Scotch is a VCF file that lists all breakpoints discovered, their confidence, and their type. 
