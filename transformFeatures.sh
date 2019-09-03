#!/bin/bash
#Transforms stdin feature matrix

meanDepth="$1"
meanNReads="$2"

#Make sure they aren't 0
#need to use bc to compare floats
if (( $( echo "$meanDepth == 0 || $meanNReads == 0" | bc -l ) ))
then
	echo Exiting: mean depth or mean nReads is 0
	exit 1
fi

#what to set feautre to in case nReads is 0 and we try to divide by it
divByZero=1000

#apply transformation
awk='BEGIN { 

#need for calculating delta features later
lastDepth = 0; 
lastNReads = 0;
lastMapQ = 0;
lastBaseQ = 0;
lastAllSC = 0; 
lastEdgeSC = 0;
lastIns = 0;
lastEdgeDel = 0; 
lastSCCons = 0;
lastSCQual = 0;

} { 

#	#$1: chrom, $2: pos -- stay unchanged
#	#$3: depth, $4: nReads -- are normalized against meanDepth and meanNReads

depthNorm = $3 / meanDepth; 
nReadsNorm  = $4 / meanNReads;

#	#$5: mapQ, $6: baseQ -- stay unchanged
#	#$7: allSC, $8: edgeSC, $9: ins, $10: edgeDel -- record as fractions against nReads
#	#$11: scCons, $12: scQual
#	#$13: scDist

##Compute delta features
depthDelta = $3 - lastDepth; 
nReadsDelta = $4 - lastNReads;
mapQDelta = $5 - lastMapQ;
baseQDelta = $6 - lastBaseQ;
allSCDelta = $7 - lastAllSC;
edgeSCDelta = $8 - lastEdgeSC;
insDelta = $9 - lastIns;
edgeDelDelta = $10 - lastEdgeDel;
scConsDelta = $11 - lastSCCons;
scQualDelta = $12 - lastSCQual;

meanNReadsPair = ($4 + lastNReads) / 2; 

#before changing current values, store  to compute delta features next time
lastDepth = $3; 
lastNReads = $4; 
lastMapQ = $5;
lastBaseQ = $6;
lastAllSC = $7;
lastEdgeSC = $8;
lastIns = $9;
lastEdgeDel = $10;
lastSCCons = $11;
lastSCQual = $12;

##Normalize allSC, edgeSC, ins, edgeDel against nReads
#check not dividing by zero
if ($4 != 0) { 
	$7 /= $4;
	$8 /= $4;
	$9 /= $4; 
	$10 /= $4; 
} else { 
	$7 = divByZero;
	$8 = divByZero; 
	$9 = divByZero;
	$10 = divByZero; 
}

if (meanNReadsPair != 0) {
	depthDelta /= meanNReadsPair;
	nReadsDelta /= meanNReadsPair;
	mapQDelta /= meanNReadsPair;
	baseQDelta /= meanNReadsPair;
	allSCDelta /= meanNReadsPair;
	edgeSCDelta /= meanNReadsPair;
	insDelta /= meanNReadsPair;
	edgeDelDelta /= meanNReadsPair;
} else { 
	depthDelta = divByZero;
	nReadsDelta = divByZero;
	mapQDelta = divByZero;
	baseQDelta = divByZero;
	allSCDelta = divByZero;
	edgeSCDelta = divByZero;
	insDelta = divByZero;
	edgeDelDelta = divByZero;

	#don"t need to normalize scConsDelta since it"s already [0,1]
	#dont"t need to normalize scQualDelta either
	#and really don"t need to normalize mapQDelta or baseQDelta either!!
}

##Print tab-delimited
OFS="\t"; 
print $1,$2,depthNorm,nReadsNorm,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,depthDelta,nReadsDelta,mapQDelta,baseQDelta,allSCDelta,edgeSCDelta,insDelta,edgeDelDelta,scConsDelta,scQualDelta;

}'

awk -F'\t' -v meanDepth="$meanDepth" -v meanNReads="$meanNReads" -v divByZero="$divByZero" "$awk"
