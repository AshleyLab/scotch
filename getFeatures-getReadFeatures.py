#!/usr/bin/python
#Processes reads, parsing CIGAR strings and extracting qualities
#and writes gz output
import csv
import gzip
import os
import numpy
import pysam
import re
import sys

bam = sys.argv[1]
bed = sys.argv[2]
outFeature = sys.argv[3]

#Open BAM file with pysam
bamData = pysam.AlignmentFile(bam, "rb")

#Open BED file and read regions to interrogate
regions = []
with open(bed, "r") as b:
	for region in csv.reader(b, delimiter="\t"):
		#add 1 since BED is 0-based
		regions.append([region[0], int(region[1]) + 1, int(region[2]) + 1]) 

#Prepare to write TSV GZ
outFeatureRaw = outFeature + ".raw" #when finish, rename to remove .raw
o = gzip.open(outFeatureRaw, "wt", newline="")
writer = csv.writer(o, delimiter="\t")

#where we'll store info about six features
#two qualities
#	mapQ - mapping quality
#	baseQ - base quality
#some from CIGAR string
#	allSC - base is SC
#	edgeSC - SC ends on this base
#	ins - insertion
#	edgeDel - deletion
#	scCons - SC consistency
#	meanSCQual - baseQ of SC bases
#	nHQual - proportion of reads with mapQ above threshold (13) [really calc'ed from baseQ]
#	scDist - signed mean distance along reads to SC (0 if base is SC)
#keys are [chrom][pos][featureType]
features = {}

def recordFeature(chrom, pos, value, featureType): 
	
	#value could be None (for increment allSC, edgeSC, ins, edgeDel), 
	#a float (for qualities)
	#or list [scPos, scBase] for scCons

	if chrom not in features: 
		features[chrom] = {}

	if pos not in features[chrom]: 
		features[chrom][pos] = {}
	
	if featureType == "scCons": 
		
		scPos = value[0]
		scBase = value[1]

		#features[chrom][pos][featureType] is all for scCons
		#a dict, where keys are positions
			#and values are dictionaries
				#of keys A, T, C, G with values how many times each occurs

		if featureType not in features[chrom][pos]: 
			features[chrom][pos][featureType] = {}

		if scPos not in features[chrom][pos][featureType]: 
			features[chrom][pos][featureType][scPos] = {}

		if scBase not in features[chrom][pos][featureType][scPos]:
			features[chrom][pos][featureType][scPos][scBase] = 1
		else: 
			features[chrom][pos][featureType][scPos][scBase] += 1

	elif featureType in ["allSC", "edgeSC", "ins", "edgeDel"]:

		#set to 1 or increment
		if featureType not in features[chrom][pos]: 
			features[chrom][pos][featureType] = 1
		else: 
			features[chrom][pos][featureType] += 1

	else: #mapQ, baseQ, scQual, nHQual, scDist
		#append to array so can take mean later
		if featureType not in features[chrom][pos]: 
			features[chrom][pos][featureType] = [value]
		else: 
			features[chrom][pos][featureType].append(value)

def recordInsertion(i, chrom, oStart, oEnd, _):
	
	#mark the insertion point
	recordFeature(chrom, oStart, None, "ins")

def recordDeletion(i, chrom, oStart, oEnd, _):

	#mark the edges of the deletion
	recordFeature(chrom, oStart, None, "edgeDel")

	if oEnd != oStart - 1: 
		#the deletion is not a single base
		#(in which case there's only position to mark)
		#actually, for single-base deletions, should we mark the base twice? 

		#decrement to avoid off by 1
		recordFeature(chrom, oEnd - 1, None, "edgeDel")

def queryPosForEdgeSCPos(scOnLeft, edgeSCPos, read):

	#get_aligned_pairs() returns list of [queryPos, pos] 
	#	except pos is None for SC bases
	#	so find queryPos for base just after/before SC
	aligned = read.get_aligned_pairs()

	#subtract one to match genomic coordinate to pysam coordinate
	boundaryPos = edgeSCPos - 1

	if scOnLeft: 
		#add 1 since need first non SC base
		boundaryPos += 1
	else: 
		#subtract 1 since need last non SC base
		boundaryPos -= 1
	
	try:
		boundaryQueryPos = [q for q in aligned if q[1] == boundaryPos][0][0]
	except:
		print("Warning: Could not find SC query position for read: {}".format(read.tostring()))
		return None

	#shift back to find query pos for edgSCPos
	if scOnLeft:
		edgeSCQueryPos = boundaryQueryPos - 1
	else: 
		edgeSCQueryPos = boundaryQueryPos + 1

	return edgeSCQueryPos

def recordSoftClipping(i, chrom, oStart, oEnd, read):

	#1. Record all SC
	#SC covers [oStart, oEnd)
	for pos in range(oStart, oEnd):
		recordFeature(chrom, pos, None, "allSC")
	
	#2. Record distance to SC
	#for positions that are not SC, record signed distance to SC
	#for positions that are SC, record 0
	rStart = read.reference_start + 1
	rEnd = rStart + read.query_length
	for pos in range(rStart, rEnd):
		if pos < oStart:
			#before SC
			distance = oStart - pos
		elif pos >= oEnd:
			#after SC
			distance = oEnd - pos
		else: 
			distance = 0
		recordFeature(chrom, pos, distance, "scDist")

	#3. Record inner edge of SC 
	scOnLeft = i == 0

	if scOnLeft:
		#last SC base is on right
		edgeSCPos = oEnd - 1
	else: 
		#last SC base is on left
		edgeSCPos = oStart

	recordFeature(chrom, edgeSCPos, None, "edgeSC")

	#3. Record base for all SC
	def recordSCBase(queryPos, pos):
		base = read.query_sequence[queryPos]
		recordFeature(chrom, edgeSCPos, [pos, base], "scCons")
		
		#also record base quality
		baseQ = read.query_qualities[queryPos]
		recordFeature(chrom, edgeSCPos, baseQ, "scQual")
	
	edgeSCQueryPos = queryPosForEdgeSCPos(scOnLeft, edgeSCPos, read)
	queryPos = edgeSCQueryPos
	if not queryPos:
		return

	if scOnLeft:
		#start at inner edge of SC, move along SC toward left of edge read
		for pos in reversed(range(oStart, oEnd)):
			recordSCBase(queryPos, pos)
			queryPos -= 1

	else: 
		#start at inner edge of SC, move along SC toward right edge of read
		for pos in range(oStart, oEnd):
			recordSCBase(queryPos, pos)
			queryPos += 1

def processCIGAR(chrom, read):

	#add 1 to make position match what's in BAM
	rStart = read.reference_start + 1
	cigar = read.cigarstring

	if not read.cigartuples: 
		return

	#loop over operations in CIGAR string
	for i, (oCode, oWidth) in enumerate(read.cigartuples):

		#we only care about 
		#	insertions (1), 
		#	deletions (2), 
		#	and soft-clipping (4)
		if oCode not in [1, 2, 4]: 
			continue

		#find the start and end of the operation
		if (i == 0): 
		
			#first operation in the read:
			#read starts at first mapped base
			#so subtract oWidth to figure out actual start
			oStart = rStart - oWidth
			oEnd = rStart

		else:
		
			#start at the start of the read,
			oStart = rStart

			#add the sum of the oWidths of 
			#all consuming operations before this
			nonConsumingCodes = [1, 4, 5, 6]

			previous = read.cigartuples[:i]
			previousConsuming = [c for c in previous if c[0] not in nonConsumingCodes]
			previousConsumingWidthsSum = sum([int(c[1]) for c in previousConsuming])

			#then adjust oStart and oEnd
			oStart += previousConsumingWidthsSum
			oEnd = oStart + oWidth

		#record operation
		{
			1 : recordInsertion, 
			2 : recordDeletion, 
			4 : recordSoftClipping

		}[oCode](i, chrom, oStart, oEnd, read)
		#pass read just for recordSoftClipping so can record SC bases

def processQualities(chrom, read):

	#separate cigar and quality extraction into separate functions? 
	#Extract qualities from reads
	#mapQ
	mapQ = read.mapping_quality
	
	#apply to all places where read matches to ref
	for pos in read.get_reference_positions():
		recordFeature(chrom, pos, mapQ, "mapQ")

	#extract base quality for each base
	for [queryPos, pos] in read.get_aligned_pairs():
		if queryPos is not None:
			#skip del, refskip?
			baseQ = read.query_qualities[queryPos]
			recordFeature(chrom, pos, baseQ, "baseQ")

def getSCCons(f, pos, nullValue):

	if "scCons" not in f or not f["scCons"]: 
		return nullValue

	#SHORTCUT: if edgeSC == 1, scCons == 1
	if "edgeSC" in f and f["edgeSC"] == 1: 
		return 1.0

	bases = f["scCons"]
	consistencies = []
	for scPos in bases: 
		scBaseCounts = bases[scPos].values()
		#this is dict, where keys are A, T, C, G
		#and values are how often each occurs
		#do values to get just occurences

		consistency = max(scBaseCounts) / sum(scBaseCounts) #ratio of times most common element appears
		consistencies.append(consistency)
	
	return numpy.mean(consistencies)

def getFeatures(chrom, pos, nullValue):

	nFeatures = 10
	if chrom not in features or pos not in features[chrom]: 
		data = [nullValue] * nFeatures
	else: 
		f = features[chrom][pos]
		scCons = getSCCons(f, pos, nullValue)

		#get meanBaseQ -- AND nHQual
		if "baseQ" not in f or not f["baseQ"]:
			meanBaseQ = nullValue
			nHQual = nullValue
		else:
			bQs = f["baseQ"]
			meanBaseQ = numpy.mean(bQs)
			nHQual = len(list(filter(lambda bQ : bQ > 13, bQs))) / float(len(bQs))
	
		#get mean mapQ
		if "mapQ" not in f or not f["mapQ"]:
			meanMapQ = nullValue
		else:
 			meanMapQ = numpy.mean(f["mapQ"])

		#get mean scQual
		if "scQual" not in f or not f["scQual"]:
			meanSCQual = nullValue
		else: 
			meanSCQual = numpy.mean(f["scQual"])

		#get mean scDist
		if "scDist" not in f or not f["scDist"]:
			meanSCDist = 1000 
			#nullValue is low, which would provide strong signal that close to SC
			#so when we don't know, use a high placeholder flag
		else: 
			meanSCDist = numpy.mean(f["scDist"])

		data = [
			meanMapQ,
			meanBaseQ,
			f["allSC"] if "allSC" in f else nullValue,
			f["edgeSC"] if "edgeSC" in f else nullValue,
			f["ins"] if "ins" in f else nullValue,
			f["edgeDel"] if "edgeDel" in f else nullValue,
			scCons,
			meanSCQual,
			nHQual,
			meanSCDist
		]

	return data
	
#Iterate over BAM file
#Populating features with values to write out
windowSize = 10000 #If bottleneck is I/O, using bigger windowSize should make faster
windowBuffer = 300
for [chrom, start, end] in regions:

	windowStart = start

	while windowStart < end: 
		
		windowEnd = windowStart + windowSize
		if windowEnd > end: 
			windowEnd = end
	
		features.clear()

		#Calculating
		#Extract reads from this region 
		#Store the results in features
		for read in bamData.fetch(chrom, max(0, windowStart - windowBuffer), windowEnd + windowBuffer):
			processCIGAR(chrom, read)
			processQualities(chrom, read)

		#Writing
		for pos in range(windowStart, windowEnd): 

			#note that this makes sure that all the positions we pull data for are in fact in this window
			#just iterating over the positions in features[chrom] could result in 
			#pulling data for positoins outside this window
			#since reads that are both inside and outside the window can be pulled by fetch

			nullValue = 0
			featuresData = getFeatures(chrom, pos, nullValue)

			writer.writerow([chrom, pos] + featuresData)

		#Shift window
		windowStart += windowSize

os.rename(outFeatureRaw, outFeature)
print("Done.")
