#!/bin/bash

# Input Parameters
input=$1 
output=$2
logfile=$3

# Run SAMConveryer Plugin

#echo "run_pipeline.pl -fork1 -SAMConverterPlugin  -i ${input} -o ${output} -endPlugin -runfork1"
run_pipeline.pl -fork1 -SAMConverterPlugin  -i ${input} -o ${output} -endPlugin -runfork1 
 
mv ${output}.log ${logfile}