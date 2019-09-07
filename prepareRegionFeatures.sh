#!/bin/bash
#Given a directory of bed files and region feature files,
#extracts those portions of the region features files described by the bed files

beds_dir="$1"
all_rfs_dir="$2"
output_trim_rfs_dir="$3"

for chrom in {1..22} X Y
do

	bed="${beds_dir}/${chrom}.bed"
	all_rfs_for_chrom="${all_rfs_dir}/${chrom}.rfs.gz"
	output_trim_rfs_for_chrom="${output_trim_rfs_dir}/${chrom}.rfs.trim.gz"

	echo "Chrom ${chrom}: ${bed} and ${all_rfs_for_chrom} --> ${output_trim_rfs_for_chrom}"

	# all_rfs_for_chrom is sorted, so posX is on line X:
	# build sed command like 
	# 	pos1,pos2p
	# 	pos3,pos4p
	#	...
	# to extract regions specified in bed file
	#	cut to extract position columns from bed
	#	awk to add 1 since bed is 0-based
	#	sed to add "," and "p"
	sed_script=$(mktemp)
	cut -f2-3 "$bed" | awk -F'\t' '{OFS=FS; $1++}1' | sed -e "s|\t|,|" -e "s|$|p|" > "$sed_script"

	# assumes $bed and $all_rfs_for_chrom have lines only from $chrom

	zcat "$all_rfs_for_chrom" | sed -n -f "$sed_script" | gzip > "$output_trim_rfs_for_chrom"

done
