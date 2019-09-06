#!/usr/bin/env python3

import csv
import gzip
import sys
from pathlib import Path
import typing

features_dir = sys.argv[1]
bed = sys.argv[2]
tmp_dir = sys.argv[3]
out_matrix = sys.argv[4]
rfs_file_path = sys.argv[5]

def get_mean_depth() -> float:
	# TODO
	return 0.0

def get_mean_nReads() -> float:
	# TODO
	return 0.0

# TODO
# actually, assert these are *not* 0, otherwise can't normalize
meanDepth: float = get_mean_depth()
meanNReads: float = get_mean_nReads()

depth_feature_gz: Path = Path(features_dir) / "depth.feat.gz"
nReads_feature_gz: Path = Path(features_dir) / "nReads.feat.gz"
read_features_gz: Path = Path(features_dir) / "read.feats.gz"
rfs_file: Path = Path(rfs_file_path)

n = 10
counter = 0

# indices
CHROM_INDEX = 0
POS_INDEX = 1

# feature indices
DEPTH_NREADS_INDEX = 2

# (in read feats)
READ_MAPQ_INDEX = 2
READ_BASEQ_INDEX = 3
READ_ALLSC_INDEX = 4
READ_EDGESC_INDEX = 5
READ_INS_INDEX = 6
READ_EDGEDEL_INDEX = 7
READ_SCCONS_INDEX = 8
READ_SCQUAL_INDEX = 9
READ_NHQUAL_INDEX = 10
READ_SCDIST_INDEX = 11

# for each position, keep track of previous's features to compute deltas
lastDepth, lastNReads, lastMapQ, lastBaseQ, lastAllSC, lastEdgeSC, lastIns, lastEdgeDel, lastSCCons, lastSCQual = (0, 0, 0, 0, 0, 0, 0, 0, 0, 0)

# what to set normalized feature to, in case nReads is 0 at a spot
DIV_BY_ZERO = 1000

with gzip.open(depth_feature_gz, "rt") as depth_f, \
	gzip.open(nReads_feature_gz, "rt") as nReads_f, \
	gzip.open(read_features_gz, "rt") as read_f: #, \
	#open(rfs_file) as rfs_f:

	for depth_line, nReads_line, read_line in \
		zip(csv.reader(depth_f, delimiter="\t"), csv.reader(nReads_f, delimiter="\t"), csv.reader(read_f, delimiter="\t")): 
		print(depth_line, nReads_line, read_line)

		# assert position1 == position2 == position3 == position4
		chrom: str = str(depth_line[CHROM_INDEX])
		assert chrom == nReads_line[CHROM_INDEX] == read_line[CHROM_INDEX], f"Chroms don't match: {depth_line[CHROM_INDEX]}, {nReads_line[CHROM_INDEX]}, {read_line[CHROM_INDEX]}"
		pos: int = int(depth_line[CHROM_INDEX])
		assert depth_line[POS_INDEX] == nReads_line[POS_INDEX] == read_line[POS_INDEX], f"Positions don't match: {depth_line[POS_INDEX]}, {nReads_line[POS_INDEX]}, {read_line[POS_INDEX]}"

		# combine
		depth: float = float(depth_line[DEPTH_NREADS_INDEX])
		nReads: float = float(nReads_line[DEPTH_NREADS_INDEX])
		read_feats: list = read_line[READ_MAPQ_INDEX:(READ_SCDIST_INDEX + 1)]
		combined = [chrom, pos, depth, nReads] + read_feats

		# transform
		# chrom, pos -- stay unchanged
		# depth, nReads are normalized against meanDepth and meanNReas
		depthNorm: float = float(depth) / meanDepth
		nReaadsNorm: float = float(depth) / meanNReads

		# compute delta features
		depthDelta = depth - lastDepth
		nReadsDelta = nReads - lastNReads
		mapQDelta = read_feats[READ_MAPQ_INDEX] - lastMapQ
		baseQDelta = read_feats[READ_BASEQ_INDEX] - lastBaseQ
		allSCDelta = read_feats[READ_ALLSC_INDEX] - lastAllSC
		edgeSCDelta = read_feats[READ_EDGESC_INDEX] - lastEdgeSC
		insDelta = read_feats[READ_INS_INDEX] - lastIns
		edgeDelDelta = read_feats[READ_EDGEDEL_INDEX] - lastEdgeDel
		scConsDelta = read_feats[READ_SCCONS_INDEX] - lastSCCons
		scQualDelta = read_feats[READ_SCQUAL_INDEX] - lastSCQual
		
		# nHQual, scDist -- stay unchanged

		# capture before update lastNReads
		meanNReadsPair = (nReads + lastNReads) / 2
	
		# update last values
		lastDepth = depth
		lastNReads = nReads
		lastMapQ = read_feats[READ_MAPQ_INDEX]
		lastBaseQ = read_feats[READ_BASEQ_INDEX]
		lastAllSC = read_feats[READ_ALLSC_INDEX]
		lastEdgeSC = read_feats[READ_EDGESC_INDEX]
		lastIns = read_feats[READ_INS_INDEX]
		lastEdgeDel = read_feats[READ_EDGEDEL_INDEX]
		lastSCCons = read_feats[READ_SCCONS_INDEX]
		lastSCQual = read_feats[READ_SCQUAL_INDEX

		# normalize allSC, edgeSC, ins, edgeDel against nReads
		# TODO: make sure nReads is float
		if nReads != 0:
			allSC = read_feats[READ_ALLSC_INDEX] / nReads
			edgeSC = read_feats[READ_EDGESC_INDEX] / nReads
			ins = read_feats[READ_INS_INDEX] / nReads
			edgeDel = read_feats[READ_EDGEDEL_INDEX] / nReads
		else:
			allSC = edgeSC = ins = edgeDel = DIV_BY_ZERO

		# for delta features, normalize against meanNReads
		if meanNReadsPair != 0:
			mNRP: float = float(meanNReadsPair)
			depthDelta /= mNRP
			nReadsDelta /= mNRP
			mapQDelta /= mNRP # TODO: don't need to normalize this or baseQ?
			baseQDelta /= mNRP
			allSCDelta /= mNRP
			edgeSCDelta /= mNRP
			insDelta /= mNRP
			edgeDelDelta /= mNRP
		else:
			depthDelta = nReadsDelta = mapQDelta = baseQDelta = allSCDelta = edgeSCDelta = insDelta = edgeDelDelta = DIV_BY_ZERO

		transformed = [chrom, pos, depthNorm, nReasdNorm, 

		if counter > n:
			break
		counter += 1
