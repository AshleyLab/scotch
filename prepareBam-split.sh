#!/bin/bash

rmdupBam="$1"
splitBam="$2"
bed="$3"

samtools view -b "$rmdupBam" -L "$bed" > "$splitBam"
samtools index "$splitBam"
