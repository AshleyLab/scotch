#!/bin/bash
#Calculate mean of feature file (depth.feat or nReads.feat)

awk -F"\t" '{s += $3} END {print s / NR}' "$1"
