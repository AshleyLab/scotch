#!/bin/bash
#Combine .feat{,s} files from "$featuresDir"

featuresDir="$1"
featureSet="$2"

#read.feats should be: mapQ, baseQ, ...
paste <(zcat "${featuresDir}/depth.feat.gz") \
	<(zcat "${featuresDir}/nReads.feat.gz" | cut -f3-) \
	<(zcat "${featuresDir}/read.feats.gz" | cut -f3-)
