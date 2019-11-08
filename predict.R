#!/usr/bin/env Rscript

#Configure
options(warn = 1, scipen = 999, keep.source = TRUE, error = 
  quote({ 

    cat("Environment:\n", file=stderr()); 
    dump.frames(); #writes to last.dump

    n <- length(last.dump)
    calls <- names(last.dump)
    cat(paste("  ", 1L:n, ": ", calls, sep = ""), sep = "\n", file=stderr())
    cat("\n", file=stderr())

    if (!interactive()) {
      q()
    }
  }))

library(data.table)
library(reshape)
library(randomForest)

getHeader <- function() { 

	return(c("chrom","pos","depthNorm","nReadsNorm","mapQ","baseQ","allSCRatio","edgeSCRatio","insRatio","edgeDelRatio","scCons","scQual","nHQual","scDist","depthDeltaNorm","nReadsDeltaNorm","mapQDeltaNorm","baseQDeltaNorm","allSCDeltaNorm","edgeSCDeltaNorm","insDeltaNorm","edgeDelDeltaNorm","scConsDelta","scQualDelta","nistHC","repMasker","segDups","LCR","gc50","gc1000","map100","uniq35","depthDeltaNorm2","nReadsDeltaNorm2","mapQDeltaNorm2","baseQDeltaNorm2","allSCDeltaNorm2","edgeSCDeltaNorm2","insDeltaNorm2","edgeDelDeltaNorm2","scConsDelta2","scQualDelta2"))

}

#Save results from prediction
saveResults <- function(results, outPath) {

	write.table(format(results, digits=3), outPath, quote=F, row.names=F, col.names=F, append=shouldAppend, sep="\t")
}

#make predictions
doPredict <- function(fm, model, outPath) {

	probs = predict(model, fm[, !colnames(fm) %in% c("chrom", "pos")], type="prob")

	#get final classes too
	predClass = colnames(probs)[max.col(probs, ties.method = "first")]

	#combine probabilites (rounding to 3 places), into comma-delimited string
	allResults = data.frame(fm$chrom, fm$pos, predClass, probs) #probsString)

	#remove normal positions
	results = allResults[predClass != "n", ]
	
	saveResults(results, outPath)
}

#Read arguments
args = commandArgs(trailingOnly=TRUE)
modelPath = args[1]
outPath = args[2]

#Read in model
model = readRDS(modelPath)

windowStart = 0
windowSize = 300000
shouldAppend = FALSE

#Read in labeled feature matrix in WINDOWS
print(paste("Reading..."))
input <- file("stdin")
open(input)
while(length(line <- readLines(input, n=windowSize, warn=FALSE)) > 0) {

	fm <- setDF(colsplit(line, "\t", names=getHeader()))
	
	print(paste("Read", nrow(fm), "lines."))
	
	#run rf model
	print("Running predict...")
	doPredict(fm, model, outPath)

	windowStart = windowStart + windowSize
	shouldAppend = TRUE
}

print("Done.")
