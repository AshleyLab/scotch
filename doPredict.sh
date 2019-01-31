#!/bin/bash

fmPath="$1"
outPath="$2"
modelPath="$3"

scotchDir=$(dirname "$0")

predictR="${scotchDir}/predict.R"
zcat "$fmPath" | Rscript "$predictR" "$modelPath" "$outPath"
