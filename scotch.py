#!/usr/bin/env python3

import argparse
from pathlib import Path
import subprocess
import sys
import typing

# Is the split-bam stage even really necessary? 
# Let's try removing it. 

# constants
CHROMS = list(str(c) for c in range(1, 23)) + ["X", "Y"]
AUX_DIR = "aux"
MODEL_NAME = "scotch-r1-wgs-model.RData"
GATK_JAR_NAME = "gatk-3.8.jar"

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

def get_tmp_dir(project_dir: Path) -> Path:
	return project_dir / "tmp/"

# bam convenience funcs
def get_rmdup_bam(project_dir: Path) -> Path:
	return get_bams_dir(project_dir) / "rmdup.bam"

def get_unclip_bam(project_dir: Path, for_chrom: str) -> Path:
	return get_bams_dir(project_dir, for_chrom) / f"{for_chrom}.unclip.bam"

# feature convenience funcs
def get_depth_feature_gz(project_dir: Path, for_chrom: str) -> Path:
	return get_features_dir(project_dir, for_chrom) / "depth.feat.gz"

def get_depth_feature_stats(project_dir: Path, for_chrom: str) -> Path:
	return get_features_dir(project_dir, for_chrom) / "depth.feat.stats"

def get_nReads_feature_gz(project_dir: Path, for_chrom: str) -> Path:
	return get_features_dir(project_dir, for_chrom) / "nReads.feat.gz"

def get_nReads_feature_stats(project_dir: Path, for_chrom: str) -> Path:
	return get_features_dir(project_dir, for_chrom) / "nReads.feat.stats"

def get_read_features_gz(project_dir: Path, for_chrom: str) -> Path:
	return get_features_dir(project_dir, for_chrom) / "read.feats.gz"

def get_feature_matrix_gz(project_dir: Path, for_chrom: str) -> Path:
	return get_features_dir(project_dir, for_chrom) / "matrix.txt.gz"

# results convenience funcs
def get_tsv_results(project_dir: Path, for_chrom: str) -> Path:
	assert for_chrom in CHROMS, f"Unrecognized chrom {for_chrom}"
	return get_results_dir(project_dir) / f"results.{for_chrom}.tsv"

def get_vcf_results_stub(project_dir: Path, for_chrom: str) -> Path:
	assert for_chrom in CHROMS, f"Unrecognized chrom {for_chrom}"
	return get_results_dir(project_dir) / f"results.{for_chrom}"

# bed convenience funcs
def get_bed_file(beds_dir: Path, for_chrom: str) -> Path:
	assert for_chrom in CHROMS, f"Unrecognized chrom {for_chrom}" 
	return beds_dir / f"{for_chrom}.bed"

# rfs convenience func
def get_all_rfs_file(all_rfs_dir: Path, for_chrom: str) -> Path:
	assert for_chrom in CHROMS, f"Unrecognized chrom {for_chrom}"
	return all_rfs_dir / f"{for_chrom}.rfs.gz"

def get_trim_rfs_file(trim_rfs_dir: Path, for_chrom: str) -> Path:
	assert for_chrom in CHROMS, f"Unrecognized chrom {for_chrom}"
	return trim_rfs_dir / f"{for_chrom}.rfs.trim.gz"

# check bam is valid
def quickcheck_bam(bam: Path) -> bool:
	bam_path: str = str(bam.resolve())
	return subprocess.call(["samtools", "quickcheck", bam_path]) == 0

# check GATK version is compatible
def check_gatk_38(gatk_jar: Path) -> bool:
	check_cmd: List = ["java", "-jar", str(gatk_jar.resolve()), "--version"]
	version_info: str = subprocess.check_output(check_cmd).decode()
	print(f"GATK version info: {version_info}")
	return version_info.startswith("3.8")

# run a script
def run_script(script_name, *args):
	scotch_dir: Path = Path(__file__).parent
	script: Path = scotch_dir / script_name
	str_args = [str(a) for a in args]
	print(f"running {script} with {str_args}")

	# TODO; make this nicer
	if script.suffix != ".py":
		subprocess.call([script] + str_args)
	else:
		subprocess.call(["python", script] + str_args)

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
	# unclip-bam
	def get_status_for_bam(bam: Path) -> str:
		if bam.is_file() and quickcheck_bam(bam):
			return DONE
		elif bam.is_file():
			return INVALID
		else:
			return INCOMPLETE
	unclip_bam_row = [get_status_for_bam(get_unclip_bam(project_dir, c)) for c in CHROMS]
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

def unclip_bam(args):
	# args we've already validated
	chrom: str = args.chrom
	project_dir: Path = Path(args.project_dir)
	bed_file: Path = get_bed_file(Path(args.beds_dir), chrom)

	# validate args specific to this stage: fasta_ref, gatk_jar, gatk_mem
	assert args.fasta_ref, "fasta_ref must be specified for unclip-bam stage"
	fasta_ref: Path = Path(args.fasta_ref)
	assert fasta_ref.is_file(), "fasta_ref must be a file that exists"

	assert args.gatk_jar, "gatk_jar must be specified for unclip-bam stage"
	gatk_jar: Path = Path(args.gatk_jar)
	assert gatk_jar.is_file(), "gatk_jar must be a file that exists"
	assert check_gatk_38(gatk_jar), "Scotch requires GATK 3.8"

	assert isinstance(args.gatk_mem, int), "gatk_mem must be an int"
	gatk_mem: int = int(args.gatk_mem)
	
	# prepare tmp dir
	tmp_dir: Path = get_tmp_dir(project_dir)
	tmp_dir.mkdir(exist_ok=True)

	# validate results of last stage: rmdup-bam
	rmdup_bam: Path = get_rmdup_bam(project_dir)
	assert rmdup_bam.is_file(), f"unclip-bam --chrom={chrom} needs to access {rmdup_bam}: try running rmdup-bam"
	assert quickcheck_bam(rmdup_bam), f"{rmdup_bam} is invalid: try running rmdup-bam"

	# get path to new bam
	bams_dir_for_chrom = get_bams_dir(project_dir, chrom)
	bams_dir_for_chrom.mkdir(exist_ok=True)
	unclip_bam: Path = get_unclip_bam(project_dir, chrom)
	assert not unclip_bam.exists(), f"unclip-bam --chrom={chrom} writes to {unclip_bam} but that already exists, please delete or stash"

	# run pipeline script
	script_name: str = "prepareBam-unclip.sh"
	run_script(script_name, rmdup_bam, bed_file, fasta_ref, gatk_jar, tmp_dir, gatk_mem, unclip_bam)

def get_features_depth(args):
	# args we've already validated
	chrom: str = args.chrom
	project_dir: Path = Path(args.project_dir)
	bed_file: Path = get_bed_file(Path(args.beds_dir), chrom)

	# validate args specific to this stage: fasta_ref, gatk_jar, gatk_mem
	assert args.fasta_ref, "fasta_ref must be specified for unclip-bam stage"
	fasta_ref: Path = Path(args.fasta_ref)
	assert fasta_ref.is_file(), "fasta_ref must be a file that exists"

	assert args.gatk_jar, "gatk_jar must be specified for unclip-bam stage"
	gatk_jar: Path = Path(args.gatk_jar)
	assert gatk_jar.is_file(), "gatk_jar must be a file that exists"
	assert check_gatk_38(gatk_jar), "Scotch requires GATK 3.8"

	assert isinstance(args.gatk_mem, int), "gatk_mem must be an int"
	gatk_mem: int = int(args.gatk_mem)
	
	# prepare tmp dir
	tmp_dir: Path = get_tmp_dir(project_dir)
	tmp_dir.mkdir(exist_ok=True)

	# validate results of last stage: rmdup-bam
	rmdup_bam: Path = get_rmdup_bam(project_dir)
	assert rmdup_bam.is_file(), f"unclip-bam --chrom={chrom} needs to access {rmdup_bam}: try running rmdup-bam"
	assert quickcheck_bam(rmdup_bam), f"{rmdup_bam} is invalid: try running rmdup-bam"

	# get path to new feature
	get_features_dir(project_dir).mkdir(exist_ok=True)
	get_features_dir(project_dir, chrom).mkdir(exist_ok=True)
	depth_feature_gz: Path = get_depth_feature_gz(project_dir, chrom)
	assert not depth_feature_gz.exists(), f"get-features-depth --chrom={chrom} writes to {depth_feature_gz} but that already exists, please delete or move"

	# path to feature stats
	depth_feature_stats: Path = get_depth_feature_stats(project_dir, chrom)	

	# run pipeline script
	script_name: str = "getFeatures-getReadCount.sh"
	run_script(script_name, rmdup_bam, bed_file, fasta_ref, gatk_jar, tmp_dir, gatk_mem, depth_feature_gz, depth_feature_stats)

def get_features_nReads(args):
	# args we've already validated
	chrom: str = args.chrom
	project_dir: Path = Path(args.project_dir)
	bed_file: Path = get_bed_file(Path(args.beds_dir), chrom)

	# validate args specific to this stage: fasta_ref, gatk_jar, gatk_mem
	assert args.fasta_ref, "fasta_ref must be specified for unclip-bam stage"
	fasta_ref: Path = Path(args.fasta_ref)
	assert fasta_ref.is_file(), "fasta_ref must be a file that exists"

	assert args.gatk_jar, "gatk_jar must be specified for unclip-bam stage"
	gatk_jar: Path = Path(args.gatk_jar)
	assert gatk_jar.is_file(), "gatk_jar must be a file that exists"
	assert check_gatk_38(gatk_jar), "Scotch requires GATK 3.8"

	assert isinstance(args.gatk_mem, int), "gatk_mem must be an int"
	gatk_mem: int = int(args.gatk_mem)
	
	# prepare tmp dir
	tmp_dir: Path = get_tmp_dir(project_dir)
	tmp_dir.mkdir(exist_ok=True)

	# validate restuls of last stage: unclip-bam
	unclip_bam: Path = get_unclip_bam(project_dir, chrom)
	assert unclip_bam.is_file(), f"get-features-nReads --chrom={chrom} needs to access {unclip_bam}: try running unclip-bam --chrom={chrom}"
	assert quickcheck_bam(unclip_bam), f"{unclip_bam} is invalid: try running unclip-bam --chrom={chrom}"
	
	# get path to new feature
	get_features_dir(project_dir).mkdir(exist_ok=True)
	get_features_dir(project_dir, chrom).mkdir(exist_ok=True)
	nReads_feature_gz: Path = get_nReads_feature_gz(project_dir, chrom)
	assert not nReads_feature_gz.exists(), f"get-features-nReads --chrom={chrom} writes to {nReads_feature_gz} but that already exists, please delete or move"

	# path to feature stats
	nReads_feature_stats: Path = get_nReads_feature_stats(project_dir, chrom)

	# run pipeline script
	script_name: str = "getFeatures-getReadCount.sh"
	run_script(script_name, unclip_bam, bed_file, fasta_ref, gatk_jar, tmp_dir, gatk_mem, nReads_feature_gz, nReads_feature_stats)

def get_features_read(args):
	# args we've already validated
	chrom: str = args.chrom
	project_dir: Path = Path(args.project_dir)
	bed_file: Path = get_bed_file(Path(args.beds_dir), chrom)	

	# validate results of last stage: rmdup-bam
	rmdup_bam: Path = get_rmdup_bam(project_dir)
	assert rmdup_bam.is_file(), f"get-features-read --chrom={chrom} needs to access {rmdup_bam}: try running rmdup-bam"
	assert quickcheck_bam(rmdup_bam), f"{rmdup_bam} is invalid: try running rmdup-bam"

	# get path to new feature
	get_features_dir(project_dir).mkdir(exist_ok=True)
	get_features_dir(project_dir, chrom).mkdir(exist_ok=True)
	read_features: Path = get_read_features_gz(project_dir, chrom)
	assert not read_features.exists(), f"get-features-read --chrom={chrom} writes to {read_features} but that already exists, please delete or move"

	# run pipeline script
	script_name: str = "getFeatures-getReadFeatures.sh"
	run_script(script_name, rmdup_bam, bed_file, read_features)

def compile_features(args):
	# args we've already validated
	chrom: str = args.chrom
	project_dir: Path = Path(args.project_dir)

	# validate args specific to this stage: rfs_dir
	assert args.trim_rfs_dir, "compile-features required trim_rfs_dir, a directory of region feature files (can be produced by prepare-region-features)"
	trim_rfs_dir: Path = Path(args.trim_rfs_dir)
	assert trim_rfs_dir.is_dir(), "trim_rfs_dir must be a directory that exists"
	trim_rfs_file: Path = get_trim_rfs_file(trim_rfs_dir, chrom)
	assert trim_rfs_file.is_file(), f"trim_rfs_dir must contain a file for {chrom}, {trim_rfs_file}"

	# prepare tmp dir
	tmp_dir: Path = get_tmp_dir(project_dir)
	tmp_dir.mkdir(exist_ok=True)	

	# validate input from previous stage: get-features-depth
	depth_feature: Path = get_depth_feature_gz(project_dir, chrom)
	assert depth_feature.is_file(), f"compile-features --chrom={chrom} needs access to {depth_feature}: try running get-features-depth --chrom={chrom}"

	# validate input from previous stage: get-features-nReads
	nReads_feature: Path = get_nReads_feature_gz(project_dir, chrom)
	assert nReads_feature.is_file(), f"compile-features --chrom={chrom} needs access to {nReads_feature}: try running get-features-nReads --chrom={chrom}"

	# validate input from previous stage: getfeatures-read
	read_features: Path = get_read_features_gz(project_dir, chrom)
	assert read_features.is_file(), f"compile-features --chrom={chrom} needs access to {read_features}: try running get-features-read --chrom={chrom}"

	# TODO: assert exist and n loci processed match
	depth_feature_stats: Path = get_depth_feature_stats(project_dir, chrom)
	nReads_feature_stats: Path = get_nReads_feature_stats(project_dir, chrom)	

	# get path to feature matrix
	features_dir: Path = get_features_dir(project_dir, chrom)
	feature_matrix: Path = get_feature_matrix_gz(project_dir, chrom)
	assert not feature_matrix.exists(), f"compile-features --chrom={chrom} writes to {feature_matrix} but that already exists, please delete or move"

	# run pipeline script
	script_name: str = "compileFeatures.py"
	run_script(script_name, depth_feature, nReads_feature, read_features, trim_rfs_file, depth_feature_stats, nReads_feature_stats, feature_matrix)

def predict(args):
	# args we've already validated
	chrom: str = args.chrom
	project_dir: Path = Path(args.project_dir)

	model_path: Path = Path(__file__).parent / AUX_DIR / MODEL_NAME

	# validate args specific to this stage: model_path, fasta_ref
	assert args.fasta_ref, "fasta_ref must be specified for predict"
	fasta_ref: Path = Path(args.fasta_ref)
	assert fasta_ref.is_file(), "fasta_ref must be a file that exists"

	# validate input from previous stage: compile-features
	feature_matrix: Path = get_feature_matrix_gz(project_dir, chrom)
	assert feature_matrix.is_file(), "predict --chrom=chrom needs access to {feature_matrix}: try running compile-features --chrom={chrom}"

	# get output files for this stage
	get_results_dir(project_dir).mkdir(exist_ok=True)
	tsv_results: Path = get_tsv_results(project_dir, chrom)
	assert not tsv_results.exists(), f"predict --chrom={chrom} writes to {tsv_results} but that already exists, please delete or move"
	vcf_results_stub: Path = get_vcf_results_stub(project_dir, chrom)
	# assert actual children do not exist ?
	#assert not vcf_results_stub.exists(), f"predict --chrom={chrom} writes to {vcf_results} but that already exists, please delete or move"

	# run pipeline script
	script_name: str = "doPredict.sh"
	run_script(script_name, feature_matrix, tsv_results, model_path, fasta_ref, vcf_results_stub)

def prepare_region_features(args):
	# args we've already validated
	beds_dir: Path = Path(args.beds_dir)
	
	# validate args specific to this stage: beds_dir, all_rfs_dir, output_trim_rfs_dir
	assert args.all_rfs_dir, "prepare-region-features requires all_rfs_dir (available at https://github.com/AshleyLab/scotch-data)"
	all_rfs_dir: Path = Path(args.all_rfs_dir)
	assert all_rfs_dir.is_dir(), "all_rfs_dir must be a directory that exists"
	
	assert args.output_trim_rfs_dir, "prepare-region-features requires output_trim_rfs_dir (can be an empty directory) to know where to write output"
	output_trim_rfs_dir: Path = Path(args.output_trim_rfs_dir)
	assert output_trim_rfs_dir.is_dir(), "output_trim_rfs_dir must be a directory that exists"
		
	for chrom in CHROMS:
		bed_file: Path = get_bed_file(beds_dir, chrom)
		assert bed_file.is_file(), f"beds_dir must contain a file for {chrom}, {bed_file}"
		all_rfs_file = get_all_rfs_file(all_rfs_dir, chrom)
		assert all_rfs_file.is_file(), f"all_rfs_dir must contain a file for {chrom}, {all_rfs_file}"
		output_trim_rfs_file = get_trim_rfs_file(output_trim_rfs_dir, chrom)
		assert not output_trim_rfs_file.exists(), f"prepare-region-features writes to {output_trim_rfs_file} but that already exists, please delete or move"

	# run pipeline script
	script_name: str= "prepareRegionFeatures.sh"
	run_script(script_name, beds_dir, all_rfs_dir, output_trim_rfs_dir)


# map keywords to functions that execute pipeline stages
COMMANDS = {
	"rmdup-bam" : rmdup_bam,
	"unclip-bam": unclip_bam,
	"get-features-depth": get_features_depth,
	"get-features-nReads": get_features_nReads,
	"get-features-read": get_features_read,
	"compile-features": compile_features,
	"predict": predict,
	"prepare-region-features": prepare_region_features,
	"check-status": check_status
}

if __name__ == "__main__":

	# get GATK JAR in aux/ directory
	aux_gatk_jar_path: Path = Path(__file__).parent / AUX_DIR / GATK_JAR_NAME
	aux_gatk_jar: str = str(aux_gatk_jar_path.resolve())

	parser = argparse.ArgumentParser(description="Process args")
	parser.add_argument("command", choices=COMMANDS.keys(), type=str, help="Stage of Scotch pipeline to execute")
	parser.add_argument("--project_dir", type=str, help="Absolute path to directory where Scotch files are stored")
	parser.add_argument("--bam", type=str, help="BAM on which to execute pipeline, required for initial stage, rmdup-bam")
	parser.add_argument("--chrom", choices=CHROMS, type=str, help="Chrom to execute stage of pipeline on, required for all stages of pipeline except rmdup-bam and prepare-region-features")
	parser.add_argument("--beds_dir", type=str, help="Absolute path to directory where one bed file for each chrom [1-22, X, Y] is stored as ${chrom}.bed")
	parser.add_argument("--fasta_ref", type=str, help="Absolute path to genome build FASTA reference, required for unclip-bam, get-features-depth, get-features-nReads, compile-features")
	parser.add_argument("--gatk_jar", type=str, default=aux_gatk_jar, help="Absolute path to GATK 3.8 JAR file, required for unclip-bam, get-features-depth, get-features-nReads, default in aux/")
	parser.add_argument("--gatk_mem", type=int, default=5, help="Amount of memory, in GB, to run GATK with (default: 5)")
	parser.add_argument("--all_rfs_dir", type=str, help="For prepare-region-features: absolute path to directory with region features for all positions(available at https://github.com/AshleyLab/scotch-data)")
	parser.add_argument("--output_trim_rfs_dir", type=str, help="For prepare-region-features: absolute path to directory to output trimmed region features to")
	parser.add_argument("--trim_rfs_dir", type=str, help="Absolute path to directory with one region feature file for each chrom [1-22, X, Y], obtained from prepare-region-features")

	args = parser.parse_args()

	command: str = args.command
	if command not in ["prepare-region-features"]:
		# validate project_dir
		assert args.project_dir, "All pipeline stages except prepare-region-features require project_dir"
		project_dir: Path = Path(args.project_dir)
		assert project_dir.is_dir(), "project_dir must be a directory that exists"

	if command not in ["rmdup-bam", "check-status", "prepare-region-features"]:
		# validate chrom	
		chrom: str = args.chrom
		assert chrom in CHROMS, "chrom [1-22, X, Y] must be specified for stages other than rmdup-bam"

	if command not in ["rmdup-bam", "predict", "check-status"]:
		# validate beds_dir
		assert args.beds_dir, "beds_dir must be specified for stages other than rmdup-bam and predict"
		beds_dir: Path = Path(args.beds_dir)
		assert beds_dir.is_dir(), "beds_dir must be a directory that exists"
		if args.chrom:
			bed_file: Path = get_bed_file(beds_dir, chrom)
			assert bed_file.is_file(), "beds_dir must contain a bed file for the specified chrom" 
		
	COMMANDS[command](args)
