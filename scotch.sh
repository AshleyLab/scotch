#!/bin/bash
#Prepare BAM for Scotch
#Remove PCR duplicates, split by chromosome

workingDir="$1"
bedsDir="$2"
bam="$3"
id="$4"
fastaRef="$5"
gatkJAR="$6"
rfsDir="$7"

#Constants
mem="4"

#directory with Scotch scripts
scotchDir=$(dirname "$0")
#train-22-dOne-N/models/model.epsilon.N1.iii.3.RData
modelPath="${scotchDir}/wgs-model.RData" 

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
rmdupScript="${scotchDir}/prepareBam-rmdup.sh"

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
	splitScript="${scotchDir}/prepareBam-split.sh"

	script="${doDir}/split.${chrom}.sh"	
	echo \#!/bin/bash > "$script"
	echo sh "$splitScript" "$rmdupBam" "$splitBam" "$bed" >> "$script"

	#B. Make a version of each split bam where soft-clipped bases are normal
	#This must happen AFTER $splitBam is created in split.${chrom}.sh
	unclipSplitBam="${bamsDir}/${chrom}.unclip.bam"
	unclipScript="${scotchDir}/prepareBam-unclip.sh"

	script="${doDir}/unclip.${chrom}.sh"
	echo \#!/bin/bash > "$script"
	echo sh "$unclipScript" "$splitBam" "$bed" "$fastaRef" "$gatkJAR" "$tmpDir" "$mem" "$unclipSplitBam" >> "$script"

	##2. Make Features
	#A. Depth
	getReadCount="${scotchDir}/getFeatures-getReadCount.sh"
	outFeature="${featuresDir}/depth.feat"

	script="${doDir}/gd.${chrom}.sh"
	echo \#!/bin/bash > "$script"
	echo sh "$getReadCount" "$splitBam" "$bed" "$fastaRef" "$gatkJAR" "$tmpDir" "$mem" "$outFeature" >> "$script"

	#B. nReads
	getReadCount="${scotchDir}/getFeatures-getReadCount.sh"
	outFeature="${featuresDir}/nReads.feat"
	
	script="${doDir}/gn.${chrom}.sh"
	echo \#!/bin/bash > "$script"
	echo sh "$getReadCount" "$unclipSplitBam" "$bed" "$fastaRef" "$gatkJAR" "$tmpDir" "$mem" "$outFeature" >> "$script"

	#C. Read Features
	getReadFeatures="${scotchDir}/getFeatures-getReadFeatures.sh"
	outFeature="${featuresDir}/read.feats"
		
	script="${doDir}/gr.${chrom}.sh"
	echo \#!/bin/bash > "$script"
	echo sh "$getReadFeatures" "$splitBam" "$bed" "$outFeature" >> "$script"
		
	##3. Compile Features
	compileFeatures="${scotchDir}/compileFeatures.sh"
	outMatrix="${featuresDir}/matrix.${chrom}.txt"
	
	script="${doDir}/compile.${chrom}.sh"
	echo \#!/bin/bash > "$script" 
	echo bash "$compileFeatures" "$featuresDir" "$bed" "$chrom" "$tmpDir" "$outMatrix" "$rfsDir" >> "$script"

	##4. Run Random Forest
	doPredict="${scotchDir}/doPredict.sh"
	outResults="${resultsDir}/results.${chrom}.scotch.tsv"
	vcfResults="${resultsDir}/results.${chrom}.vcf"

	script="${doDir}/predict.${chrom}.sh"
	echo \#!/bin/bash > "$script" 
	echo bash "$doPredict" "${outMatrix}.gz" "$outResults" "$modelPath" "$fastaRef" "$vcfResults" >> "$script" 
done
