#!/usr/bin/python

import csv, os, pysam
variantsPath = sys.argv[1]
fastaPath = sys.argv[2]
outPath = sys.argv[3]

#print header
def printHeader(): 
	h = "##fileformat=VCFv4.1"
	"""
	##fileformat=VCFv4.1
	##phasing=none
    ##INDIVIDUAL=TRUTH
    ##SAMPLE=<ID=TRUTH,Individual="TRUTH",Description="bamsurgeon spike-in">
    ##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
	##ALT=<ID=DEL_L, Description="Deletion Start"
	##ALT=<ID=DEL_R, Description="Deletion Right">
    ##ALT=<ID=INS,Description="Insertion">
    ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
    #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSPIKEIN""")

#process variants
def processVariant(v):

	[chrom, pos, predType, probs] = v
	reference = fasta.fetch(chrom, int(pos) - 1, int(pos)).upper()

	if predType == "ins":
		alt = "<INS>"
	elif predType == "del_L": 
		alt = "<DEL_L>"
	elif predType == "del_R": 
		alt = "<DEL_R>"
	elif predType == "dOne": 
	
		#base at this pos is one deleted, shift back one
		alternate = reference	
		reference = fasta.fetch(chrom, int(pos) - 2, int(pos)).upper()
		pos--
	else: 
		print("ERROR: Variant has unexpected type {}".format(predType))


	writer.writerow([chrom, start, ".", reference, alternate]) #more?

#Execution starts here
#Read in FASTA
fasta = pysam.FastaFile(fastaPath)

#Read in variants
with open(variantsPath, "r") as v: 
	for variant in csv.reader(v, delimiter="\t"):
		processVariant(v)
