#!/bin/bash

# Input parameter
input=$1
min_tag_number=$2
fastq=$3  

# Run TagCountToFastq
run_pipeline.pl -fork1 -TagCountToFastqPlugin  -i ${input} -o ${fastq} -c ${min_tag_number} -endPlugin -runfork1 
