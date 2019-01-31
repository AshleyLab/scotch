#!/bin/bash
#Removes PCR-duplicates

bam="$1"
rmdupBam="$2"

samtools rmdup "$bam" "$rmdupBam"
samtools index "$rmdupBam"
