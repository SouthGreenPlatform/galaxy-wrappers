#!/bin/bash
input_exonerate=$1
database=$2
input_region=$3
output_exonerate=$4
/usr/local/bioinfo/galaxy/galaxy_dist/tools/sequence_comparisons/exonerate.pl $input_exonerate $database $input_region;
mv "$1.exonerate" $4;

