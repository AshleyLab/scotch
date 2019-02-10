#!/usr/bin/python

import csv, os, pysam, sys
variantsPath = sys.argv[1]
fastaPath = sys.argv[2]

#print header
def writeHeader():

	h = ["##fileformat=VCFv4.1"]
	h.append("##phasing=none")
	h.append("##ALT=<ID=DEL_L, Description=\"Deletion Start\">")
	h.append("##ALT=<ID=DEL_R, Description=\"Deletion Right\">")
	h.append("##ALT=<ID=INS,Description=\"Insertion\">")
	h.append("##INFO=<ID=PROBS,Number=1,Type=String,Description=\"Class probabilites from random forest model\">")
	h.append("##FORMAT=<ID=GT,Number=1, Type=String,Description=\"Genotype\">")
	h.append("\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLE"]))

	for headerItem in h: 
		print(headerItem)

#process variants
def processVariant(v):

	[chrom, p, predType, prob1, prob2, prob3, prob4, prob5] = v
	pos = int(p)

	probs = ",".join([prob1, prob2, prob3, prob4, prob5])
	ref = fasta.fetch(chrom, pos - 1, pos).upper()

	if predType == "ins":
		alt = "<INS>"
	elif predType == "del_L": 
		alt = "<DEL_L>"
	elif predType == "del_R": 
		alt = "<DEL_R>"
	elif predType == "dOne": 
	
		#base at this pos is one deleted, shift back one
		alt = ref
		ref = fasta.fetch(chrom, int(pos) - 2, int(pos)).upper()
		pos -= 1
	else: 
		sys.stderr.write("ERROR: Variant has unexpected type {}".format(predType))
		return
	
	ID     = "."
	QUAL   = "100"
	FILTER = "PASS"
	INFO   = "PROBS=" + probs
	FORMAT = "GT"
	SAMPLE = "./."
	#      "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLE"
	data = [chrom, pos, ID, ref, alt, QUAL, FILTER, INFO, FORMAT, SAMPLE]
	print("\t".join(str(x) for x in data))

#Execution starts here
#Read in FASTA
fasta = pysam.FastaFile(fastaPath)

writeHeader()

#Read in variants
with open(variantsPath, "r") as v: 
	for variant in csv.reader(v, delimiter="\t"):
		processVariant(variant)
