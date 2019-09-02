#!/bin/bash
#Calculates depth of coverage (includes or excludes SC bases 
#depending on whether called with unclipped bam or not from getFeatures.sh)
#and writes, gzipped, to outFeature

bam="$1"
bed="$2"
fastaRef="$3"
gatkJAR="$4"
tmpDir="$5"
mem="$6"
outFeature="$7"

#depthDoc="${outFeature}.doc"
log="${outFeature}.log"

#Use GATK to calculate depth across portion of BAM in BED
java -Djava.io.tmpdir="$tmpDir" -Xmx"$mem"g -jar "$gatkJAR" \
	-T DepthOfCoverage \
	-mbq 0 -mmq 0 \
	--omitIntervalStatistics \
	--omitLocusTable \
	--omitPerSampleStats \
	--includeRefNSites \
	-L "$bed" \
	-I "$bam" \
	-o /dev/stdout \
	-R "$fastaRef" 2> "$log" | \
grep -vi warn	`#0` | \
tail -n +2	`#1` | \
cut -f-2 	`#2` | \
sed "s|:|\t|"   `#3` | \
gzip > "$outFeature"

#       1. Use tail to remove header
#       2. Use cut to extract chrom, pos, value columns
#       3. Use sed to parse into chrom\tpos\value format
#added --includeRefNSites to make sure emits output for ambiguous bases (because read.feats will have)
echo Done.
