#!/bin/bash

import csv
import gzip
from pathlib import Path
import sys
from typing import Any, List, NamedTuple

# constants
# flag value for when nReads is 0 and we try to normalize against it
DIV_BY_ZERO = 1000
# zero-indexed column in .stats file for feat mean value
STATS_MEAN_INDEX: int = 1
# delimiter expected in .feat and .stats files
DELIMITER = "\t"

class Coordinates(NamedTuple):
	""" Packages chrom, pos pair to represent genomic coordinates """
	chrom: str
	pos: int

class RegionFeats(NamedTuple):
	""" Collection of eight region features that describe for a genomic locus 
	    stored as strings because not modified, used in any calculations,
	    and want to avoid rounding/padding
	"""
	nistHC: str
	repMasker: str
	segDups: str
	LCR: str
	gc50: str
	gc1000: str
	map100: str
	uniq35: str

class BaseFeats(NamedTuple):
	""" Collection of base features, input directly to this script """
	coords: Coordinates
	depth: float  # from depth.feat
	nReads: float  # from nReads.feat
	mapQ: float  # from read.feats ...
	baseQ: float
	allSC: float
	edgeSC: float
	ins: float
	edgeDel: float
	scCons: float
	scQual: float
	nHQual: float
	scDist: float
	region_feats: RegionFeats # from .rfs file

class NormalizedBaseFeats(NamedTuple):
	""" Collection of normalized features calculated from BaseFeats """
	coords: Coordinates
	normalized_depth: float
	normalized_nReads: float
	mapQ: float
	baseQ: float
	normalized_allSC: float
	normalized_edgeSC: float
	normalized_ins: float
	normalized_edgeDel: float
	scCons: float
	scQual: float
	nHQual: float
	scDist: float
	region_feats: RegionFeats

class DeltaFeats(NamedTuple):
	""" Collection of "delta" features that describe the changes between two BaseFeats """
	coords: Coordinates
	last_coords: Coordinates
	mean_nReads: float # used to normalize
	depth_delta: float
	nReads_delta: float
	mapQ_delta: float
	baseQ_delta: float
	allSC_delta: float
	edgeSC_delta: float
	ins_delta: float
	edgeDel_delta: float
	scCons_delta: float
	scQual_delta: float

class NormalizedDeltaFeats(NamedTuple): 
	""" Collection of normalized "delta" features calculated from DeltaFeats """
	coords: Coordinates
	last_coords: Coordinates
	normalized_depth_delta: float
	normalized_nReads_delta: float
	normalized_mapQ_delta: float
	normalized_baseQ_delta: float
	normalized_allSC_delta: float
	normalized_edgeSC_delta: float
	normalized_ins_delta: float
	normalized_edgeDel_delta: float
	scCons_delta: float
	scQual_delta: float

def get_mean_from_stats(stats_file: Path) -> float: 
	with open(stats_file) as stats:
		lines = stats.readlines()
	return float(lines[0].strip().split(DELIMITER)[STATS_MEAN_INDEX])

# confirm chrom and pos match across feature files, and return
def get_coordinates(depth_line: List, nReads_line: List, read_line: List, rfs_line: List) -> Coordinates:
	
	assert depth_line[0] == nReads_line[0] == read_line[0] == rfs_line[0], "Chromosomes do not match: {depth_line[0]}, {nReads_line[0]}, {read_line[0]}, {rfs_line[0]}"
	assert depth_line[1] == nReads_line[1] == read_line[1] == rfs_line[1], "Positions do not match: {depth_line[1]}, {nReads_line[1]}, {read_line[1]}, {rfs_line[1]}"
	return Coordinates(chrom=depth_line[0], pos=int(depth_line[1]))

# compute base feats from feature files
def get_base_feats(coords: Coordinates, depth_line: List, nReads_line: List, read_line: List, rfs_line: List) -> BaseFeats:

	# need str, int wrapper?	
	return BaseFeats(
		coords=coords,
		depth=float(depth_line[2]),
		nReads=float(nReads_line[2]),
		mapQ=float(read_line[2]),
		baseQ=float(read_line[3]),
		allSC=float(read_line[4]),
		edgeSC=float(read_line[5]),
		ins=float(read_line[6]),
		edgeDel=float(read_line[7]),
		scCons=float(read_line[8]),
		scQual=float(read_line[9]),
		nHQual=float(read_line[10]),
		scDist=float(read_line[11]),
		region_feats=RegionFeats(
			nistHC=rfs_line[2],
			repMasker=rfs_line[3],
			segDups=rfs_line[4],
			LCR=rfs_line[5],
			gc50=rfs_line[6],
			gc1000=rfs_line[7],
			map100=rfs_line[8],
			uniq35=rfs_line[9]
		)
	)

# compute normalized feats from BaseFeats
def get_normalized_feats(base_feats: BaseFeats, matrix_mean_depth: float, matrix_mean_nReads: float) -> NormalizedBaseFeats:

	# normalize depth, nReads against matrix-wide means
	normalized_depth: float = base_feats.depth / matrix_mean_depth
	normalized_nReads: float = base_feats.nReads / matrix_mean_nReads

	# normalize allSC, edgeSC, ins, edgeDel against nReads	
	nReads: float = base_feats.nReads
	def normalize(feat: float) -> float:
		if nReads != 0:
			return feat / nReads
		return DIV_BY_ZERO

	normalized_allSC = normalize(base_feats.allSC)
	normalized_edgeSC = normalize(base_feats.edgeSC)
	normalized_ins = normalize(base_feats.ins)
	normalized_edgeDel = normalize(base_feats.edgeDel)

	return NormalizedBaseFeats(
		coords=base_feats.coords,
		normalized_depth=normalized_depth,
		normalized_nReads=normalized_nReads,
		mapQ=base_feats.mapQ,
		baseQ=base_feats.baseQ,
		normalized_allSC=normalized_allSC,
		normalized_edgeSC=normalized_edgeSC,
		normalized_ins=normalized_ins,
		normalized_edgeDel=normalized_edgeDel,
		scCons=base_feats.scCons,
		scQual=base_feats.scQual,
		nHQual=base_feats.nHQual,
		scDist=base_feats.scDist,
		region_feats=base_feats.region_feats
	)	

# compute changes in features between two BaseFeats
def get_delta_feats(base_feats: BaseFeats, last_base_feats: BaseFeats) -> DeltaFeats:
	
	# used in get_normalized_delta_feats
	mean_nReads: float = (base_feats.nReads + last_base_feats.nReads) / 2

	return DeltaFeats(
		coords=base_feats.coords,
		last_coords=last_base_feats.coords,
		mean_nReads=mean_nReads,
		depth_delta=base_feats.depth - last_base_feats.depth,
		nReads_delta=base_feats.nReads - last_base_feats.nReads,
		mapQ_delta=base_feats.mapQ - last_base_feats.mapQ,
		baseQ_delta=base_feats.baseQ - last_base_feats.baseQ,
		allSC_delta=base_feats.allSC - last_base_feats.allSC,
		edgeSC_delta=base_feats.edgeSC - last_base_feats.edgeSC,
		ins_delta=base_feats.ins - last_base_feats.ins,
		edgeDel_delta=base_feats.edgeDel - last_base_feats.edgeDel,
		scCons_delta=base_feats.scCons - last_base_feats.scCons,
		scQual_delta=base_feats.scQual - last_base_feats.scQual
		
	)

# compute normalized feats from DeltaFeats
def get_normalized_delta_feats(delta_feats: DeltaFeats) -> NormalizedDeltaFeats:

	# normalize depth, nReads_, mapQ_, baseQ_, allSC_, edgeSC_, ins_, edgeDel_delta
	# against mean nReads from two loci this tracks the delta between
	mean_nReads: float =  delta_feats.mean_nReads
	def normalize(feat: float) -> float:
		if mean_nReads != 0:
			return feat / mean_nReads
		return DIV_BY_ZERO

	normalized_depth_delta = normalize(delta_feats.depth_delta)
	normalized_nReads_delta = normalize(delta_feats.nReads_delta)
	normalized_mapQ_delta = normalize(delta_feats.mapQ_delta)
	normalized_baseQ_delta = normalize(delta_feats.baseQ_delta)
	normalized_allSC_delta = normalize(delta_feats.allSC_delta)
	normalized_edgeSC_delta = normalize(delta_feats.edgeSC_delta)
	normalized_ins_delta = normalize(delta_feats.ins_delta)
	normalized_edgeDel_delta = normalize(delta_feats.edgeDel_delta)

	return NormalizedDeltaFeats(
		coords=delta_feats.coords,
		last_coords=delta_feats.last_coords,
		normalized_depth_delta=normalized_depth_delta,
		normalized_nReads_delta=normalized_nReads_delta,
		normalized_mapQ_delta=normalized_mapQ_delta,
		normalized_baseQ_delta=normalized_baseQ_delta,
		normalized_allSC_delta=normalized_allSC_delta,
		normalized_edgeSC_delta=normalized_edgeSC_delta,
		normalized_ins_delta=normalized_ins_delta,
		normalized_edgeDel_delta=normalized_edgeDel_delta,
		scCons_delta=delta_feats.scCons_delta,
		scQual_delta=delta_feats.scQual_delta
	)

# combine features into a single list
def collate_feats(last_normalized_base_feats: NormalizedBaseFeats, 
		last_normalized_delta_feats: NormalizedDeltaFeats, 
		normalized_delta_feats: NormalizedDeltaFeats) -> List:

	last_nbf: NormalizedBaseFeats = last_normalized_base_feats
	last_ndf: NormalizedDeltaFeats = last_normalized_delta_feats
	ndf: NormalizeDeltaFeats = normalized_delta_feats

	return [last_nbf.coords.chrom, last_nbf.coords.pos,
		last_nbf.normalized_depth, last_nbf.normalized_nReads,
		last_nbf.mapQ, last_nbf.baseQ,
		last_nbf.normalized_allSC, last_nbf.normalized_edgeSC, last_nbf.normalized_ins,last_nbf.normalized_edgeDel,
		last_nbf.scCons, last_nbf.scQual, last_nbf.nHQual, last_nbf.scDist,
		last_ndf.normalized_depth_delta, last_ndf.normalized_nReads_delta, last_ndf.normalized_mapQ_delta, last_ndf.normalized_baseQ_delta,
		last_ndf.normalized_allSC_delta, last_ndf.normalized_edgeSC_delta, last_ndf.normalized_ins_delta, last_ndf.normalized_edgeDel_delta,
		last_ndf.scCons_delta, last_ndf.scQual_delta,
		last_nbf.region_feats.nistHC, last_nbf.region_feats.repMasker, last_nbf.region_feats.segDups, last_nbf.region_feats.LCR,
		last_nbf.region_feats.gc50, last_nbf.region_feats.gc1000, last_nbf.region_feats.map100, last_nbf.region_feats.uniq35,
		ndf.normalized_depth_delta, ndf.normalized_nReads_delta, ndf.normalized_mapQ_delta, ndf.normalized_baseQ_delta,
		ndf.normalized_allSC_delta, ndf.normalized_edgeSC_delta, ndf.normalized_ins_delta, ndf.normalized_edgeDel_delta,
		ndf.scCons_delta, ndf.scQual_delta]

def format_feat(feat) -> str:
	if isinstance(feat, float):
		return f"{feat:0.6g}"
	return feat

# write feats through csv writer
def write_feats(writer: Any, feats: List) -> None:
	writer.writerow(format_feat(feat) for feat in feats)

# get csv readers for reading csv files
def get_reader(feature: Any) -> Any:
	return csv.reader(feature, delimiter=DELIMITER)

# null versions of classes we defined, used as placeholders
EMPTY_COORDS = Coordinates(chrom="0", pos=0)
EMPTY_REGION_FEATS = RegionFeats(nistHC=0.0, repMasker=0.0, segDups=0.0, LCR=0.0, gc50=0.0, gc1000=0.0, map100=0.0, uniq35=0.0)
EMPTY_BASE_FEATS = BaseFeats(coords=EMPTY_COORDS, depth=0.0, nReads=0.0, mapQ=0.0, baseQ=0.0, allSC=0.0, edgeSC=0.0, ins=0.0, edgeDel=0.0,
	scCons=0.0, scQual=0.0, nHQual=0.0, scDist=0.0, region_feats=EMPTY_REGION_FEATS)
EMPTY_NORMALIZED_BASE_FEATS = NormalizedBaseFeats(coords=EMPTY_COORDS, normalized_depth=0.0, normalized_nReads=0.0, mapQ=0.0, baseQ=0.0, 
	normalized_allSC=0.0,normalized_edgeSC=0.0,normalized_ins=0.0,normalized_edgeDel=0.0,scCons=0.0,scQual=0.0,nHQual=0.0,scDist=0.0,region_feats=EMPTY_REGION_FEATS)
EMPTY_DELTA_FEATS = DeltaFeats(coords=EMPTY_COORDS, last_coords=EMPTY_COORDS, mean_nReads=0.0, depth_delta=0.0, nReads_delta=0.0, mapQ_delta=0.0, baseQ_delta=0.0,
	allSC_delta=0.0, edgeSC_delta=0.0, ins_delta=0.0, edgeDel_delta=0.0, scCons_delta=0.0, scQual_delta=0.0)
EMPTY_NORMALIZED_DELTA_FEATS = NormalizedDeltaFeats(coords=EMPTY_COORDS, last_coords=EMPTY_COORDS, normalized_depth_delta=0.0, normalized_nReads_delta=0.0, 
	normalized_mapQ_delta=0.0, normalized_baseQ_delta=0.0, normalized_allSC_delta=0.0, normalized_edgeSC_delta=0.0,
	normalized_ins_delta=0.0, normalized_edgeDel_delta=0.0, scCons_delta=0.0, scQual_delta=0.0)

counter = 0
increment = 10000
increment_hit = 0

if __name__ == "__main__":

	# parse args
	depth_feature_gz_path = sys.argv[1]
	nReads_feature_gz_path = sys.argv[2]
	read_features_gz_path = sys.argv[3]
	region_features_path = sys.argv[4]
	depth_feature_stats = sys.argv[5]
	nReads_feature_stats = sys.argv[6]
	out_matrix_gz_path = sys.argv[7]

	# set up output
	out_matrix_gz = gzip.open(out_matrix_gz_path, "wt", newline="")
	out_writer = csv.writer(out_matrix_gz, delimiter="\t", lineterminator="\n")

	# get mean depth, nReads per sample from .stats files
	# decide whether to use path, make it nicer
	matrix_mean_depth: float = get_mean_from_stats(Path(depth_feature_stats))
	matrix_mean_nReads: float = get_mean_from_stats(Path(nReads_feature_stats))
	assert matrix_mean_depth != 0 and matrix_mean_nReads != 0, "Mean sample depth or nReads must not be 0"
	print(f"matrix_mean_depth is {matrix_mean_depth} and matrix_mean_nReads is {matrix_mean_nReads}")

	# iterate over feature files
	with gzip.open(depth_feature_gz_path, "rt") as depth_feature, \
		gzip.open(nReads_feature_gz_path, "rt") as nReads_feature, \
		gzip.open(read_features_gz_path, "rt") as read_features, \
		open(region_features_path) as region_features:

		is_first_line = True
		features_iterator = zip(get_reader(depth_feature), get_reader(nReads_feature),
			get_reader(read_features), get_reader(region_features))

		# store features of previous row
		last_base_feats: BaseFeats = EMPTY_BASE_FEATS
		last_normalized_base_feats: NormalizedBaseFeats = EMPTY_NORMALIZED_BASE_FEATS
		last_normalized_delta_feats: NormalizedDeltaFeats = EMPTY_NORMALIZED_DELTA_FEATS

		# iterate over feature files
		for depth_line, nReads_line, read_line, rfs_line in features_iterator:
			
			#print(f"depth_line: {depth_line}, nReads_line: {nReads_line}, read_line: {read_line}, rfs_line: {rfs_line}")	

			# get base features from this row
			coords: Coordinates = get_coordinates(depth_line, nReads_line, read_line, rfs_line)
			base_feats: BaseFeats = get_base_feats(coords, depth_line, nReads_line, read_line, rfs_line)
			normalized_base_feats: NormalizedBaseFeats = get_normalized_feats(base_feats, matrix_mean_depth, matrix_mean_nReads)

			# compute delta features against previous row
			delta_feats: DeltaFeats = get_delta_feats(base_feats, last_base_feats)
			normalized_delta_feats: NormalizedDeltaFeats = get_normalized_delta_feats(delta_feats)

			# write features for previous row, with these delta features and its own
			if not is_first_line:
				last_collated_feats: List = collate_feats(last_normalized_base_feats, last_normalized_delta_feats, normalized_delta_feats)
				write_feats(out_writer, last_collated_feats)

			# save this row's features for next pass
			is_first_line = False
			last_base_feats = base_feats
			last_normalized_base_feats = normalized_base_feats
			last_normalized_delta_feats = normalized_delta_feats

			counter += 1
			if counter > increment_hit:
				print(f"Processed {counter} lines.")
				increment_hit += increment

		# collate, write one more time
		collated_feats: List = collate_feats(last_normalized_base_feats, last_normalized_delta_feats, EMPTY_NORMALIZED_DELTA_FEATS)
		write_feats(out_writer, collated_feats)

	print("Done.")



