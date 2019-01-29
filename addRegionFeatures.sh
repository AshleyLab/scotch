#!/bin/bash
#Accepts feature matrix through stdin
#And adds region features

#use to subset region features
#assume has lines from one chrom only
bed="$1"

if [[ "$bed" == "EMPTY" ]]
then
	cat -
fi

#list of 8 region features for all positions in chrom
rfs="$2"

tmpDir="$3"
#$rfs is sorted, so posX is on line X:
#build sed command like 
#pos1,pos2p;pos3,pos4p;pos5,pos6p;...
#	cut to extract position columns from bed
#	awk to add 1 since  bed is 0-based
#	sed to replace \t --> ,
#	tr to replace newlines with ; (tr can only replace one character)
#	sed to replace ; --> p;
sedScript=$(mktemp -p "$tmpDir")
cut -f2-3 "$bed" | awk -F'\t' '{OFS=FS; $1++}1' | sed -e "s|\t|,|" -e "s|$|p|" > "$sedScript"

#Paste in region features
#	sed to subset appropriate lines from $rfs
#	cut to remove chrom and pos from $rfs
paste - <(sed -n -f "$sedScript" "$rfs" | cut -f3-)
