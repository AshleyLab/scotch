#!/bin/bash

featuresDir="$1"
bed="$2"
chrom="$3"
tmpDir="$4"
outMatrix="$5"
rfsDir="$6"

#directory with Scotch scripts
scotchDir=$(dirname "$0")

#1. Combine
combine="${scotchDir}/combineFeatures.sh"

#2. Transform
#calculate means for normalizing
meanDepth=$("${scotchDir}"/calcMean.sh <(zcat "$featuresDir"/depth.feat.gz))
meanNReads=$("${scotchDir}"/calcMean.sh <(zcat "$featuresDir"/nReads.feat.gz))
transform="${scotchDir}/transformFeatures.sh"

#3. Add RFs
add="${scotchDir}/addRegionFeatures.sh"
rfs="${rfsDir}/${chrom}.rfs"

bash "$combine" "$featuresDir" "$featureSet" | "$transform" "$meanDepth" "$meanNReads" | "$add" "$bed" "$rfs" "$tmpDir" > "${outMatrix}.no_slide"

#4. Slide
slide="${scotchDir}/slideFeatures.sh"

bash "$slide" "${outMatrix}.no_slide" "${outMatrix}"
rm "${outMatrix}.no_slide"
