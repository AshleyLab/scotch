#!/bin/bash
#c

bam="$1"
bed="$2"
outFeature="$3"

#script to call
s="${SCOTCH}/scripts/final/getFeatures-getReadFeatures.py"
python "$s" "$bam" "$bed" "$outFeature"
