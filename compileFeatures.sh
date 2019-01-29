#!/bin/bash

featuresDir="$1"
bed="$2"
chrom="$3"
tmpDir="$4"
outMatrix="$5"

#1. Combine
combine="$SCOTCH/scripts/final/combineFeatures.sh"

#2. Transform
#calculate means for normalizing
meanDepth=$("$SCOTCH"/scripts/core/calcMean.sh <(zcat "$featuresDir"/depth.feat.gz))
meanNReads=$("$SCOTCH"/scripts/core/calcMean.sh <(zcat "$featuresDir"/nReads.feat.gz))
transform="$SCOTCH/scripts/final/transformFeatures.sh"

#3. Add RFs
add="$SCOTCH/scripts/final/addRegionFeatures.sh"
rfs="${SCOTCH}/aux/region-features-labeled-shift/${chrom}.rfs"

bash "$combine" "$featuresDir" "$featureSet" | "$transform" "$meanDepth" "$meanNReads" | "$add" "$bed" "$rfs" "$tmpDir" > "${outMatrix}.no_slide"

#4. Slide
slide="${SCOTCH}/scripts/final/slideFeatures.sh"

bash "$slide" "${outMatrix}.no_slide" "${outMatrix}"
rm "${outMatrix}.no_slide"
