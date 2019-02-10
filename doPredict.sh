#!/bin/bash

fmPath="$1"
outPath="$2"
modelPath="$3"
fastaRef="$4"
vcfResults="$5"

scotchDir=$(dirname "$0")

#run random forest
predictR="${scotchDir}/predict.R"
zcat "$fmPath" | Rscript "$predictR" "$modelPath" "$outPath"

#convert results to VCF
makeVCF="${scotchDir}/makeVCF.py"
python "$makeVCF" "$outPath" "$fastaRef" > "$vcfResults"
