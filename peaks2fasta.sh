#!/bin/bash
# JRA 2022

# Check the number of command line arguments
if [ $# -ne 2 ]; then
        script_name=$(basename $0)
        echo "Usage: $script_name peaks.bed reference.fasta"
        exit 1
fi

INPUT=$(basename $1)
REFERENCE=$2

# remove header, keep first 3 columns, sort by coordinate
grep -v "#" $INPUT | cut -f1-3 | sort -k1,1 -k2,2n > $INPUT"_tmp"

# pull out the underlying sequences
bedtools getfasta -fi $REFERENCE -bed $INPUT"_tmp" -fo "${INPUT%.*}".fasta
rm $INPUT"_tmp"
