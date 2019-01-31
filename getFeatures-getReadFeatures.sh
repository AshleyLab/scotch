#!/bin/bash

bam="$1"
bed="$2"
outFeature="$3"

#script to call
s="${scotchDir}/getFeatures-getReadFeatures.py"
python "$s" "$bam" "$bed" "$outFeature"

gzip "$outFeature"
