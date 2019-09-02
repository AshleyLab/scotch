#!/usr/bin/env python3

import argparse
from pathlib import Path
import subprocess
import sys
import typing

# constants
CHROMS = list(str(c) for c in range(1, 23)) + ["X", "Y"]

# paths used by multiple stages
def get_bams_dir(project_dir: Path, for_chrom: str = None) -> Path:
	if for_chrom:
		assert for_chrom in CHROMS, f"Unrecognized chrom {for_chrom}" 
		return project_dir / "bams/" / for_chrom
	else:
		return project_dir / "bams/"

def get_features_dir(project_dir: Path, for_chrom: str = None) -> Path:
	if for_chrom:
		assert for_chrom in CHROMS, f"Unrecognized chrom {for_chrom}" 
		return project_dir / "features/" / for_chrom
	else:
		return project_dir / "features/"

def get_results_dir(project_dir: Path) -> Path:
	return project_dir / "results/"

def get_rmdup_bam(project_dir: Path) -> Path:
	return get_bams_dir(project_dir) / "rmdup.bam"

# bam convenience funcs
def get_split_bam(project_dir: Path, for_chrom: str) -> Path:
	return get_bams_dir(project_dir, for_chrom) / f"{for_chrom}.bam"

def get_unclip_bam(project_dir: Path, for_chrom: str) -> Path:
	return get_bams_dir(project_dir, for_chrom) / f"{for_chrom}.unclip.bam"

# feature convenience funcs
def get_depth_feature_gz(project_dir: Path, for_chrom: str) -> Path:
	return get_features_dir(project_dir, for_chrom) / "depth.feat.gz"

def get_nReads_feature_gz(project_dir: Path, for_chrom: str) -> Path:
	return get_features_dir(project_dir, for_chrom) / "nReads.feat.gz"

def get_read_features_gz(project_dir: Path, for_chrom: str) -> Path:
	return get_features_dir(project_dir, for_chrom) / "read.feats.gz"

def get_feature_matrix_gz(project_dir: Path, for_chrom: str) -> Path:
	return get_features_dir(project_dir, for_chrom) / "matrix.txt.gz"

# results convenience funcs
def get_scotch_tsv_results(project_dir: Path, for_chrom: str) -> Path:
	assert for_chrom in CHROMS, f"Unrecognized chrom {for_chrom}"
	return get_results_dir(project_dir) / f"results.{for_chrom}.scotch-tsv"

def get_vcf_results(project_dir: Path, for_chrom: str) -> Path:
	assert for_chrom in CHROMS, f"Unrecognized chrom {for_chrom}"
	return get_results_dir(project_dir) / f"results.{for_chrom}.vcf"

# bed convenience funcs
def get_bed_file(beds_dir: Path, for_chrom: str) -> Path:
	assert for_chrom in CHROMS, f"Unrecognized chrom {for_chrom}" 
	return beds_dir / f"{for_chrom}.bed"
#

# check bam is valid
def quickcheck_bam(bam: Path) -> bool:
	bam_path: str = str(bam.resolve())
	return subprocess.call(["samtools", "quickcheck", bam_path]) == 0

def run_script(script_name, *args):
	scotch_dir: Path = Path(__file__).parent
	script: Path = scotch_dir / script_name
	print(f"running {script} with {args}")
	if script_name.endswith(".sh"):
		print("calling")
		subprocess.call([script] + list(args))

# prints graphical report of status
def check_status(args) -> None: 
	# looks kinda like
	# 
	#   |   |   |   |   |   |   |   |   |   | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 2 | 2 | 2 |   |   
	# 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 0 | 1 | 2 | X | Y

	SEP: str = "|"
	DONE: str = "#" # stages of the pipeline completed ostensibly correctly
	INVALID: str = "?" # stages with clear issues
	INCOMPLETE: str = " " # stages not completed
	
	MAX_LABEL_LEN: int = 8
	def print_row(row: [str], label: str = None) -> None:
		padded_items = [f" {x} " for x in row]
		joined_items = SEP.join(padded_items)

		trimmed_label: str = (label or "")[:MAX_LABEL_LEN]
		formatted_label: str = f"{trimmed_label:8}"
		print(formatted_label + joined_items)

	# not per chrom
	project_dir: Path = Path(args.project_dir)
	rmdup_bam: Path = get_rmdup_bam(project_dir)
	if rmdup_bam.is_file() and quickcheck_bam(rmdup_bam):	
		rmdup_desc = "done"
	elif rmdup_bam.is_file():
		rmdup_desc = "invalid"
	else:
		rmdup_desc: "incomplete"
	rmdup_bam_row = f"     rmdup-bam: {rmdup_desc}"
	print(rmdup_bam_row)
	
	# chrom headers
	chroms_first_row  = [(c[0] if len(c) == 2 else " ") for c in CHROMS]
	chroms_second_row = [(c[1] if len(c) == 2 else c)  for c in CHROMS]
	print_row(chroms_first_row)
	print_row(chroms_second_row)

	# steps in pipeline
	# split-bam, unclip-bam
	def get_status_for_bam(bam: Path) -> str:
		if bam.is_file() and quickcheck_bam(bam):
			return DONE
		elif bam.is_file():
			return INVALID
		else:
			return INCOMPLETE
	split_bam_row  = [get_status_for_bam( get_split_bam(project_dir, c)) for c in CHROMS]
	unclip_bam_row = [get_status_for_bam(get_unclip_bam(project_dir, c)) for c in CHROMS]
	print_row(split_bam_row, label="split-bam")
	print_row(unclip_bam_row, label="unclip-bam")

	# get-features-depth, get-features-nReads, get-features-read
	def get_status_for_feature(feature: Path) -> str:
		if feature.is_file():
			return DONE
		else:
			return INCOMPLETE
	get_features_depth_row  = [get_status_for_feature( get_depth_feature_gz(project_dir, c)) for c in CHROMS]
	get_features_nReads_row = [get_status_for_feature(get_nReads_feature_gz(project_dir, c)) for c in CHROMS]
	get_features_read_row   = [get_status_for_feature( get_read_features_gz(project_dir, c)) for c in CHROMS]

	# compile-features
	def get_status_for_matrix(matrix: Path) -> str:
		if matrix.is_file():
			return DONE
		else:
			return INCOMPLETE
	compile_feature_row = [get_status_for_matrix(get_feature_matrix_gz(project_dir, c)) for c in CHROMS]

	# predict
	

# SCOTCH PIPELINE
def rmdup_bam(args):
	# args we've already validated
	project_dir: Path = Path(args.project_dir)

	# validate args specific to this stage: bam
	assert args.bam, "bam must be specified for rmdup-bam stage"
	input_bam: Path = Path(args.bam)
	assert input_bam.is_file(), "bam must be a file that exists"

	# make bams dir
	bams_dir: Path = get_bams_dir(project_dir)
	assert not bams_dir.exists(), f"rmdup-bam writes to {bams_dir} but that already exists, please delete or choose another project directory"
	bams_dir.mkdir()
		
	# rmdup args.bam and put the result 
	rmdup_bam: Path = get_rmdup_bam(project_dir)
	assert not rmdup_bam.exists(), f"rmdup-bam writes to {rmdup_bam} but that already exists, please delete or choose another project directory"
	
	# run pipeline script
	script_name: str = "prepareBam-rmdup.sh"
	run_script(script_name, input_bam, rmdup_bam)

def split_bam(args):
	# args we've already validated
	chrom: str = args.chrom
	project_dir: Path = Path(args.project_dir)
	bed_file: Path = get_bed_file(Path(args.beds_dir), chrom)

	# validate results of last stage: rmdup-bam
	rmdup_bam: Path = get_rmdup_bam(project_dir)
	assert rmdup_bam.is_file(), f"split-bam must be run after rmdup-bam because it needs to access {rmdup_bam}"
	assert quickcheck_bam(rmdup_bam), f"samtools quickcheck says {rmdup_bam} is invalid. Did rmdup-bam finish successfully?"
	bams_dir: Path = get_bams_dir(project_dir)
	assert bams_dir.is_dir(), f"split-bam must be run after rmdup-bam because it needs to access {bams_dir}"

	# create chrom bam directory
	bams_dir_for_chrom = get_bams_dir(project_dir, chrom)
	assert not bams_dir_for_chrom.exists(), f"split-bam writes to {bams_dir_for_chrom} but that already exists, please delete or choose another project directory"
	bams_dir_for_chrom.mkdir()

	# get path to where will write split bam
	split_bam: Path = get_split_bam(project_dir, chrom)
	assert not split_bam.exists(), f"split-bam writes to {split_bam} but that already exists, please delete or choose another project directory"

	# run pipeline script
	script_name: str = "prepareBam-split.sh"
	run_script(script_name, rmdup_bam, split_bam, bed_file)

def unclip_bam(args):
	# args we've already validated
	chrom: str = args.chrom
	project_dir: Path = Path(args.project_dir)
	bed_file: Path = get_bed_file(Path(args.beds_dir), chrom)

	# validate args specific to this stage: fasta_ref, gatk_jar, tmp_dir, mem
	fasta_ref: Path = None
	gatk_jar: Path = None
	tmp_dir: Path = None
	mem: int = None
	# TODO

	# validate results of last stage: split-bam
	split_bam: Path = get_split_bam(project_dir, chrom)
	assert split_bam.is_file(), f"unclip-bam --chrom={chrom} needs to access {split_bam}: try running split-bam --chrom={chrom}"
	assert quickcheck_bam(split_bam), f"{split_bam} is invalid: try running split-bam --chrom={chrom}"

	# get path to new bam
	unclip_bam: Path = get_unclip_bam(project_dir, chrom)
	assert not unclip_bam.exists(), f"unclip-bam --chrom={chrom} writes to {unclip_bam} but that already exists, please delete or stash"

	# run pipeline script
	script_name: str = "prepareBam-unclip.sh"
	run_script(script_name, split_bam, bed_file, fasta_ref, gatk_jar, tmp_dir, mem, unclip_bam)

def get_features_depth(args):
	# args we've already validated
	chrom: str = args.chrom
	project_dir: Path = Path(args.project_dir)
	bed_file: Path = get_bed_file(Path(args.beds_dir), chrom)

	# validate args specific to this stage: fasta_ref, gatk_jar, tmp_dir, mem
	fasta_ref: Path = None
	gatk_jar: Path = None
	tmp_dir: Path = None
	mem: int = None
	# TODO

	# validate results of last stage: split-bam
	split_bam: Path = get_split_bam(project_dir, chrom)
	assert split_bam.is_file(), f"unclip-bam --chrom={chrom} needs to access {split_bam}: try running split-bam --chrom={chrom}"
	assert quickcheck_bam(split_bam), f"{split_bam} is invalid: try running split-bam --chrom={chrom}"

	# get path to new feature
	depth_feature_gz: Path = get_depth_feature_gz(project_dir, chrom)
	assert not depth_feature_gz.exists(), f"get-features-depth --chrom={chrom} writes to {depth_feature_gz} but that already exists, please delete or move"

	# run pipeline script
	script_name: str = "getFeatures-getReadCount.sh"
	run_script(script_name, split_bam, bed_file, fasta_ref, gatk_jar, tmp_dir, mem, depth_feature_gz)

def get_features_nReads(args):
	# args we've already validated
	chrom: str = args.chrom
	project_dir: Path = Path(args.project_dir)
	bed_file: Path = get_bed_file(Path(args.beds_dir), chrom)

	# validate args specific to this stage: fasta_ref, gatk_jar, tmp_dir, mem
	fasta_ref: Path = None
	gatk_jar: Path = None
	tmp_dir: Path = None
	mem: int = None
	# TODO

	# validate restuls of last stage: unclip-bam
	unclip_bam: Path = get_unclip_bam(project_dir, chrom)
	assert unclip_bam.is_file(), f"get-features-nReads --chrom={chrom} needs to access {unclip_bam}: try running unclip-bam --chrom={chrom}"
	assert quickcheck_bam(unclip_bam), f"{unclip_bam} is invalid: try running unclip-bam --chrom={chrom}"
	
	# get path to new feature
	nReads_feature_gz: Path = get_nReads_feature_gz(project_dir, chrom)
	assert not nReads_feature_gz.exists(), f"get-features-nReads --chrom={chrom} writes to {nReads_feature_gz} but that already exists, please delete or move"

	# run pipeline script
	script_name: str = "getFeatures-getReadCount.sh"
	run_script(script_name, unclip_bam, bed_file, fasta_ref, gatk_jar, tmp_dir, meme, nReads_feature_gz)

def get_features_read(args):
	# args we've already validated
	chrom: str = args.chrom
	project_dir: Path = Path(args.project_dir)
	bed_file: Path = get_bed_file(Path(args.beds_dir), chrom)	

	# validate results of last stage: split-bam
	split_bam: Path = get_split_bam(project_dir, bam)
	assert split_bam.is_file(), f"get-features-read --chrom={chrom} needs to access {split_bam}: try running split-bam --chrom={chrom}"
	assert quickcheck_bam(split_bam), f"{split_bam} is invalid: try running split-bam --chrom={chrom}"

	# get path to new feature
	read_features: Path = get_read_features(project_dir, chrom)
	assert not read_features.exists(), f"get-features-read --chrom={chrom} writes to {read_features} but that already exists, please delete or move"

	# run pipeline script
	script_name: str = "getFeatures-getReadFeatures.sh"
	run_script(script_name, split_bam, bed_file, read_features)

def compile_features(args):
	# args we've already validated
	chrom: str = args.chrom
	project_dir: Path = Path(args.project_dir)
	bed_file: Path = get_bed_file(Path(args.beds_dir), chrom)

	# validate args specific to this stage: tmp_dir, rfs_dir
	tmp_dir: Path = None
	rfs_dir: Path = None
	# TODO

	# validate input from previous stage: get-features-depth
	depth_feature: Path = get_depth_feature(project_dir, chrom)
	assert depth_feature.is_file(), f"compile-features --chrom={chrom} needs access to {depth_feature}: try running get-features-depth --chrom={chrom}"

	# validate input from previous stage: get-features-nReads
	nReads_feature: Path = get_depth_feature(project_dir, chrom)
	assert depth_feature.is_file(), f"compile-features --chrom={chrom} needs access to {nReads_features}: try running get-features-nReads --chrom={chrom}"

	# validate input from previous stage: getfeatures-read
	read_features: Path = get_read_features(project_dir, chrom)
	assert not read_features.exists(), f"compile-features --chrom={chrom} needs access to {read_features}: try running get-features-read --chrom={chrom}"

	# get path to feature matrix
	feature_matrix: Path = get_feature_matrix(project_dir, chrom)
	features_dir: Path = get_features_dir(proejct_dir, chrom)

	# run pipeline script
	script_name: str = "compileFeatures.sh"
	run_script(script_name, features_dir, bed_file, chrom, tmp_dir, feature_matrix, rfs_dir)

def predict(args):
	# args we've already validated
	chrom: str = args.chrom
	project_dir: Path = Path(args.project_dir)
	bed_file: Path = get_bed_file(Path(args.beds_dir), chrom)

	# validate args specific to this stage: model_path, fasta_ref
	model_path: Path = None
	fasta_ref: Path = None
	
	# validate input from previous stage: compile-features
	feature_matrix: Path = get_feature_matrix(project_dir)
	assert feature_matrix.is_file(), "predict --chrom=chrom needs access to {feature_matrix}: try running compile-features --chrom={chrom}"

	# get output files for this stage
	scotch_tsv_results: Path = get_scotch_tsv_results(project_dir, chrom)
	assert not scotch_tsv_results.exists(), f"predict --chrom={chrom} writes to {scotch_tsv_results} but that already exists, please delete or move"
	vcf_results: Path = get_vcf_results(project_dir, chrom)
	assert not vcf_results.exist(), f"predict --chrom={chrom} writes to {vcf_results} but that already exists, please delete or move"

	# run pipeline script
	script_name: str = "doPredict.sh"
	run_script(script_name, feature_matrix, scotch_tsv_results, vcf_results, model_path, fasta_ref, vcf_results)

# map keywords to functions that execute pipeline stages
COMMANDS = {
	"rmdup-bam" : rmdup_bam,
	"split-bam": split_bam,
	"unclip-bam": unclip_bam,
	"get-features-depth": get_features_depth,
	"get-features-nReads": get_features_nReads,
	"get-features-read": get_features_read,
	"compile-features": compile_features,
	"predict": predict,
	"check-status": check_status
}

if __name__ == "__main__":

	parser = argparse.ArgumentParser(description="Process args")
	parser.add_argument("command", choices=COMMANDS.keys(), type=str, help="Stage of Scotch pipeline to execute")
	parser.add_argument("--project_dir", required=True, type=str, help="Absolute path to directory where Scotch files are stored")
	parser.add_argument("--bam", type=str, help="BAM on which to execute pipeline, required for initial stage, rmdup-bam")
	parser.add_argument("--chrom", choices=CHROMS, type=str, help="Chrom to execute stage of pipeline on, required for all stages of pipeline except rmdup-bam")
	parser.add_argument("--beds_dir", type=str, help="Absolute path to directory where one bed file for each chrom [1-22, X, Y] is stored as ${chrom}.bed")

	args = parser.parse_args()

	# validate project_dir
	project_dir: str = args.project_dir
	assert Path(project_dir).is_dir(), "project_dir must be a path to a directory that exists"
	
	command: str = args.command
	if command not in ["rmdup-bam", "check-status"]:
		# validate chrom	
		chrom: str = args.chrom
		assert chrom in CHROMS, "chrom must be specified for stages other than rmdup-bam"

		# validate beds_dir
		beds_dir: Path = Path(args.beds_dir)
		assert beds_dir.is_dir(), "beds_dir must be specified for stages other than rmdup-bam, and must be a directory that exists"
		bed_file: Path = get_bed_file(beds_dir, chrom)
		assert Path(bed_file).is_file(), "beds_dir must contain a bed file for the specified chrom" 
		
	COMMANDS[command](args)
