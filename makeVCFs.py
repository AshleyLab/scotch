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
INFO_TEMPLATE = "PROBS={}" 
FORMAT = "GT"
GT = "./."
ENCODE_GT = "0/1"

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
def write_variant(writers: List[Any], chrom: str, pos: int, ref: str, alt: str, info: str, gt: str) -> None:
	variant_row = [chrom, pos, ID, ref, alt, QUAL, FILTER, info, FORMAT, gt]
	for writer in writers:
		writer.writerow(variant_row)

# process variant, writing to VCFs
def process_variant(variant: List[str], writer: Any, fasta: Any) -> None:

	# unpack fields
	[chrom, raw_unshifted_pos, pred_type, prob_1, prob_2, prob_3, prob_4, prob_5] = variant
	unshifted_pos: int = int(raw_unshifted_pos)
	probs: str = ",".join([prob_1, prob_2, prob_3, prob_4, prob_5])
	info: str = INFO_TEMPLATE.format(probs)

	# shift pos, if applicate for type
	assert pred_type in PRED_TYPES, f"Variant at {chrom}:{pos} has unexpected type {pred_type}"
	shifted_pos: int = shift_pos_for_pred_type(unshifted_pos, pred_type)
	shifted_pos_ref: str = get_nucs(fasta, chrom, shifted_pos)

	# write results to standard vcf
	if pred_type == "dOne": 
		# base at shifted_pos is deleted, base at (shifted_pos - 1) is retained
		ref: str = get_nucs(fasta, chrom, shifted_pos - 1, shifted_pos + 1)
		alt: str = get_nucs(fasta, chrom, shifted_pos - 1)
		write_variant(writer, chrom, shifted_pos - 1, ref, alt, info, GT)
	else:
		alt: str = f"<{pred_type.upper()}>"
		write_variant(writer, chrom, shifted_pos, shifted_pos_ref, alt, info, GT)

if __name__ == "__main__":

	# parse args
	tsv_results_path = sys.argv[1]
	fasta_path = sys.argv[2]
	vcf_results_stub = sys.argv[3]

	# read in FASTA reference
	fasta = pysam.FastaFile(fasta_path)

	# set up output
	results_vcf = open(f"{vcf_results_stub}.vcf", "w")
	writer = csv.writer(results_vcf, delimiter=OUTPUT_DELIMITER, quoting=csv.QUOTE_NONE, quotechar="")

	# write VCF headers to output files
	chrom_lengths: Dict[str, int] = get_chrom_lengths(fasta)
	write_header(writer, chrom_lengths)

	# process variants
	with open(tsv_results_path, "r") as t: 
		for variant in csv.reader(t, delimiter="\t"):
			process_variant(variant, writer, fasta)

	# close output files
	results_vcf.close()
	for vcf in [results_vcf, encoded_del_L_results_vcf, encoded_del_R_results_vcf, encoded_ins_results_vcf, encoded_all_results_vcf]:
		vcf.close()
	
	# run_script("encode.py")
