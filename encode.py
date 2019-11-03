#!/usr/bin/env python3
# Convert Scotch output to VCF format
# Called by doPredict.sh as
# 	python makeVCF.py [tsv results] [fasta ref] [vcf results stub]
# Writes the full results in VCF format to ${output-stub}.vcf,

# Shifts the positions of del_L, dOne, and ins down by 1

import csv
import os
import pysam
import textwrap
from typing import Any, Dict, List
import typing
import sys

# constants
CHROMS = list(str(c) for c in range(1, 23)) + ["X", "Y"]

# types of Scotch calls
PRED_TYPES = ["del_L", "del_R", "dOne", "ins"]
# types of calls for which we subtract one from the position, to adjust indexing
SHIFT_TYPES = ["del_L", "dOne", "ins"]

OUTPUT_DELIMITER = "\t"
# vcf fields
ID = "."
QUAL = "100"
FILTER = "PASS"
INFO = "-"
FORMAT = "GT"
GT = "0/1"

# get chromosome lengths from fasta reference for ##contig headers
def get_chrom_lengths(fasta: Any) -> Dict[str, int]:
	return {chrom: fasta.get_reference_length(chrom) for chrom in CHROMS}

# write_header to provided csv writer
def write_header(writer: Any, chrom_lengths: Dict[str, int]) -> None:

	# generic vcf headers
	headers: [str] = ["##fileformat=VCFv4.1"]
	headers.append("##phasing=none")
	headers.append("##ALT=<ID=DEL_L, Description=\"Deletion Start\">")
	headers.append("##ALT=<ID=DEL_R, Description=\"Deletion End\">")
	headers.append("##ALT=<ID=INS,Description=\"Insertion\">")
	headers.append("##INFO=<ID=PROBS,Number=1,Type=String,Description=\"Class probabilites from random forest model\">")
	headers.append("##FORMAT=<ID=GT,Number=1, Type=String,Description=\"Genotype\">")

	for header in headers:
		writer.writerow([header])

	# contig headers
	for chrom in CHROMS:
		chrom_length: int = chrom_lengths[chrom]
		writer.writerow([f"##contig=<ID={chrom},length={chrom_length}>"])

	# column headers
	writer.writerow(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLE"])

# get nucleotides from FASTA at position
def get_nucs(ref: Any, chrom: str, start: int, end: int = None) -> str:
	if not end: 
		end = start + 1
	
	# pysam is 0-based (https://pysam.readthedocs.io/en/latest/glossary.html#term-region) so subtract 1
	zero_based_start: int = start - 1
	zero_based_end: int = end - 1
	return ref.fetch(chrom, zero_based_start, zero_based_end).upper()

# shift position of calls of certain types
def shift_pos_for_pred_type(pos: int, pred_type: str) -> int:
	if pred_type in SHIFT_TYPES:
		return pos - 1
	else:
		return pos

# for encoded vcfs, pick an arbitrary distinct allele for alt
def get_alt_for_ref(ref: str) -> str:
	ref_alt_map: Dict[str, str] = {
		"A": "T",
		"T": "C", 
		"C": "G", 
		"G": "A"
	}
	if ref in ref_alt_map: 
		return ref_alt_map[ref]
	else:
		# probably an N, ambiguous base
		return "A"

# write a variant with given fields to a list of writers
def write_variant(writers: List[Any], chrom: str, pos: int) -> None:

	ref_nuc: str = get_nucs(ref, chrom, pos)
	alt_nuc: str = get_alt_for_ref(ref_nuc)
	fields = [chrom, pos, ID, ref_nuc, alt_nuc, QUAL, FILTER, INFO, FORMAT, GT]

	for writer in writers:
		writer.writerow(fields)

# process variant, writing to VCFs
def process_variant(variant: List[str], writers: Dict[str, Any], fasta: Any) -> None:

	# unpack fields
	[chrom, str_pos, rsid, full_ref, alt_field, qual, filter_field, info, format_field, sample] = variant
	pos = int(str_pos)

	extra_fields = [info, qual, filter_field, info, format_field, sample]
	
	for alt in alts.split(","):
	
		# check whether indel
		INS_TAG = "<INS>"
		DEL_L_TAG = "<DEL_L>"
		DEL_R_TAG = "<DEL_R>" 
		DEL_TAG = "<DEL>"
		indel_tags = [INS_TAG, DEL_L_TAG, DEL_R_TAG, DEL_TAG]

		if alt in indel_tags:
	
			# alt is one of indel tags
			if (alt == INS_TAG):
			
				ins_writers = writers["ins"]
				ins_pos = pos
				write_variant(ins_writers, chrom, pos=ins_pos)

			elif (alt == DEL_L_TAG):
			
				del_L_writers = writers["del_L"]
				del_L_pos = pos + 1
				write_variant(del_L_writers, chrom=chrom, pos=del_L_pos)

			elif (alt == DEL_R_TAG):
	
				del_R_writers = writers["del_R"]
				del_R_pos = pos
				write_variant(deL_R_writers, chrom-chrom, pos=del_R_pos)

			elif (alt == DEL_TAG): # Pindel deletion

				# DEL_L
				del_L_writers = writers["del_L"]
				del_L_pos = pos + 1
				write_variant(deL_L_writers, chrom=chrom, pos=del_L_pos)

				# DEL_R
				del_R_writers = writers["del_R"] 
				# get endpoint from END tag in INFO
				end_tag = info.split(";")[0]
				del_R_pos = int(end_tag.split("=")[1])
				write_variant(del_R_writers, chrom=chrom, pos=del_R_pos)
	
		else:
			
			# alt is nucleotides, not tag
			length_diff = len(alt) - len(ref)

			if (length_diff == 0): # SNP; skip
				continue

			if (length_diff > 0): # insertion

				ins_writers = writers["ins"]
				ins_pos = pos
				write_variant(ins_writers, chrom=chrom, pos=ins_pos)

			elif (length_diff < 0): # deletion
	
				# DEL_L
				del_L_writers = writers["del_L"]
				del_L_pos = pos + 1
				write_variant(del_L_writers, chrom=chrom, pos=del_L_pos)

				# DEL_R
				del_R_writers = writers["del_R"]
				deL_R_pos = pos - length_diff
				write_variant(del_R_writers, chrom=chrom, pos=del_R_pos)

# sort numerically by position
def sort_output(output_tsv, sorted_output_tsv) -> None:

	# skip header
	skip_header_cmd = f"""awk '/^#/' {output_tsv} > {sorted_output_tsv}"""
	skip_header_output = subprocess.check_output(skip_header_cmd, shell=True)

	# sort the rest numerically by position
	sort_cmd = f"""awk '!/^#/' {output_tsv} | sort -t$'\t' -k2,2n >> {sorted_output_tsv}"""
	sort_output = subprocess.check_output(sort_cmd, shell=True)
	print(f"sort output: {sort_output}")

	# rm unsorted and check error
	rm_cmd = f"""[ "$(wc -l < {output_tsv})" -eq "$(wc -l < {sorted_output_tsv})" ] && rm {output_tsv} || echo sort error"""
	rm_output = subprocess.check_output(rm_cmd, shell=True)
	print(f"rm output: {rm_output}")

if __name__ == "__main__":

	# parse args
	vcf_results = sys.argv[1]
	fasta_path = sys.argv[2]
	vcf_results_stub = sys.argv[3]

	# read in FASTA reference
	fasta = pysam.FastaFile(fasta_path)

	# set up output
	def get_unsorted_and_sorted_paths(stub):
		return (f"{stub}.unsorted.vcf", f"{stub}.vcf")

	(unsorted_del_L_path, del_L_path) = get_unsorted_and_sorted_paths(f"{vcf_results_stub}.encode_del_L")
	(unsorted_del_R_path, del_R_path) = get_unsorted_and_sorted_paths(f"{vcf_results_stub}.encode_del_R")
	(unsorted_ins_path, ins_path) = get_unsorted_and_sorted_paths(f"{vcf_results_stub}.encode_ins")
	(unsorted_all_path, all_path) = get_unsorted_and_sorted_paths(f"{vcf_results_stub}.encode_all")

	def writer_for_vcf(vcf: Any) -> Any:
		return csv.writer(vcf, delimiter=OUTPUT_DELIMITER, quoting=csv.QUOTE_NONE, quotechar="")
	
	encoded_del_L_results_vcf = open(unsorted_del_L_path, "w")
	encoded_del_R_results_vcf = open(unsorted_del_R_path, "w")
	encoded_ins_results_vcf = open(unsorted_ins_path, "w")
	encoded_all_results_vcf = open(unsorted_all_path, "w")

	all_writer = writer_for_vcf(encoded_all_results_vcf)
	writers = {
		"del_L": [writer_for_vcf(encoded_del_L_results_vcf) + all_writer],
		"del_R": [writer_for_vcf(encoded_del_R_results_vcf) + all_writer],
		"ins": [writer_for_vcf(encoded_ins_results_vcf) + all_writer]
	}
	
	# write VCF headers to output files
	chrom_lengths: Dict[str, int] = get_chrom_lengths(fasta)
	for _, writer in writers.items():
		write_header(writer, chrom_lengths)

	# process variants
	with open(vcf_results, "r") as t: 
		for variant in csv.reader(t, delimiter="\t"):
			process_variant(variant, writers, fasta)

	# close output files
	for vcf in [encoded_del_L_results_vcf, encoded_del_R_results_vcf, encoded_ins_results_vcf, encoded_all_results_vcf]:
		vcf.close()

	# sort ouput
	# (can get out of order from multiallelic records being split)
	sort_output(unsorted_del_L_path, del_L_path)
	sort_output(unsorted_del_R_path, del_R_path)
	sort_output(unsorted_ins_path, ins_path)
	sort_output(unsorted_all_path, all_path)
	
