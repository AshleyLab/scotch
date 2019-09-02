#!/bin/bash

bam="$1"
bed="$2"
outFeature="$3"

#script to call
scotchDir=$(dirname "$0")
s="${scotchDir}/getFeatures-getReadFeatures.py"
python "$s" "$bam" "$bed" "$outFeature"
