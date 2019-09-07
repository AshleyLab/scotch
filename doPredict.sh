#!/bin/bash

fmPath="$1"
outPath="$2"
modelPath="$3"
fastaRef="$4"
vcfResultsStub="$5"

scotchDir=$(dirname "$0")

#run random forest
echo Running random forest...
predictR="${scotchDir}/predict.R"
zcat "$fmPath" | Rscript "$predictR" "$modelPath" "$outPath"

#convert results to VCF
echo Converting results to VCF...
makeVCFs="${scotchDir}/makeVCFs.py"
python "$makeVCFs" "$outPath" "$fastaRef" "$vcfResultsStub"
