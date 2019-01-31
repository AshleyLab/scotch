#!/bin/bash

fmPath="$1"
outPath="$2"
modelPath="${SCOTCH}/train-22-dOne-N/models/model.epsilon.N1.iii.3.RData"

predictR="${SCOTCH}/scripts/final/predict.R"
zcat "$fmPath" | Rscript "$predictR" "$modelPath" "$outPath"
