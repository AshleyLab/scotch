#!/bin/bash

featuresDir="$1"
bed="$2"
tmpDir="$3"
outMatrix="$4"
rfs="$5"

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
bash "$combine" "$featuresDir" "$featureSet" | "$transform" "$meanDepth" "$meanNReads" | "$add" "$bed" "$rfs" "$tmpDir" > "${outMatrix}.no_slide"

#4. Slide
slide="${scotchDir}/slideFeatures.sh"

bash "$slide" "${outMatrix}.no_slide" "${outMatrix}"
rm "${outMatrix}.no_slide"
