#!/bin/bash
#Prepare BAM for Scotch
#Remove PCR duplicates, split by chromosome

workingDir="$1"
bedsDir="$2"
bam="$3"
id="$4"

#Constants (should be read from a config file?)
mem="4"
fastaRef="/share/PI/euan/apps/bcbio/genomes/Hsapiens/GRCh37/seq/GRCh37.fa"
gatkJAR="${SCOTCH}/aux/gatk/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar"

doDir="${workingDir}/do/"
mkdir "$doDir"
tmpDir="${workingDir}/tmp/"
mkdir "$tmpDir"
bamsParentDir="${workingDir}/bams/"
mkdir "$bamsParentDir"
featuresParentDir="${workingDir}/features/"
mkdir "$featuresParentDir"
resultsDir="${workingDir}/results/"
mkdir "$resultsDir"

#0. Rmdup
rmdupToDo="${doDir}/rmdup.sh"
rmdupBam="${bamsParentDir}/${id}.rmdup.bam"
rmdupScript="${SCOTCH}/final/prepareBam-rmdup.sh"

echo \#!/bin/bash > "$rmdupToDo"
echo sh "$rmdupScript" "$bam" "$rmdupBam" >> "$rmdupToDo"

for chrom in {1..22} X Y
do
	bed="${bedsDir}/${chrom}.bed"
	bamsDir="${bamsParentDir}/${chrom}/"
	mkdir "$bamsDir"
	featuresDir="${featuresParentDir}/${chrom}/"
	mkdir "$featuresDir"

	##1. BAM Preparation
	#A. Split by chromosome
	splitBam="${bamsDir}/${chrom}.bam"
	splitScript="${SCOTCH}/final/prepareBam-split.sh"

	script="${doDir}/split.${chrom}.sh"	
	echo \#!/bin/bash > "$script"
	echo sh "$splitScript" "$rmdupBam" "$splitBam" "$bed" >> "$script"

	#B. Make a version of each split bam where soft-clipped bases are normal
	#This must happen AFTER $splitBam is created in split.${chrom}.sh
	unclipSplitBam="${SCOTCH}/scripts/final/prepareBam-unclip.sh"
	unclipScript="${bamsDir}/${chrom}.unclip.bam"

	script="${doDir}/unclip.${chrom}.sh"
	echo \#!/bin/bash > "$script"
	echo sh "$unclipScript" "$splitBam" "$bed" "$fastaRef" "$gatkJAR" "$tmpDir" "$mem" "$unclipSplitBam" >> "$script"

	##2. Make Features
	#A. Depth
	getReadCount="${SCOTCH}/scripts/final/getFeatures-getReadCount.sh"
	outFeature="${featuresDir}/depth.feat"

	script="${doDir}/gd.${chrom}.sh"
	echo \#!/bin/bash > "$script"
	echo sh "$getReadCount" "$splitBam" "$bed" "$fastaRef" "$gatkJAR" "$tmpDir" "$mem" "$outFeature" >> "$script"

	#B. nReads
	getReadCount="${SCOTCH}/scripts/final/getFeatures-getReadCount.sh"
	outFeature="${featuresDir}/nReads.feat"
	
	script="${doDir}/gn.${chrom}.sh"
	echo \#!/bin/bash > "$script"
	echo sh "$getReadCount" "$unclipSplitBam" "$bed" "$fastaRef" "$gatkJAR" "$tmpDir" "$mem" "$outFeature" >> "$script"

	#C. Read Features
	getReadFeatures="${SCOTCH}/scripts/final/getFeatures-getReadFeatures.sh"
	outFeature="${featuresDir}/read.feats"
		
	script="${doDir}/gr.${chrom}.sh"
	echo \#!/bin/bash > "$script"
	echo sh "$getReadFeatures" "$splitBam" "$bed" "$outFeature" >> "$script"
		
	##3. Compile Features
	compileFeatures="${SCOTCH}/scripts/final/compileFeatures.sh"
	outMatrix="${featuresDir}/matrix.${chrom}.txt"
	
	script="${doDir}/compile.${chrom}.sh"
	echo \#!/bin/bash > "$script" 
	echo bash "$compileFeatures" "$featuresDir" "$bed" "$chrom" "$tmpDir" "$outMatrix"  >> "$script"

	##4. Run Random Forest
	doPredict="${SCOTCH}/scripts/final/doPredict.sh"
	outResults="${resultsDir}/results.${chrom}.txt"

	script="${doDir}/predict.${chrom}.sh"
	echo \#!/bin/bash > "$script" 
	echo bash "$doPredict" "${outMatrix}.gz" "$outResults" >> "$script" 
	
	##5. Convert to VCF
done
