#!/bin/bash
#Change SC-bases to matches

childBam="$1"
bed="$2"
fastaRef="$3"
gatkJAR="$4"
tmpDir="$5"
mem="$6"
childBamUnclip="$7"

#Use GATK to change soft-clipped bases to normal bases
java -Djava.io.tmpdir="$tmpDir" -Xmx"$mem"g -jar "$gatkJAR" \
	-T ClipReads \
	-CR REVERT_SOFTCLIPPED_BASES \
	-L "$bed" \
	-I "$childBam" \
	-o "$childBamUnclip" \
	-R "$fastaRef"

echo Done.
