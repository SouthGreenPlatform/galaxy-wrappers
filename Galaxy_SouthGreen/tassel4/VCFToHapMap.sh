#!/bin/bash

vcf=$1 
hapmap=$2 
#echo "run_pipeline.pl -fork1 -vcf ${vcf} -export ${hapmap} -exportType Hapmap -runfork1"
run_pipeline.pl -fork1 -vcf ${vcf} -export ${hapmap} -exportType Hapmap -runfork1
mv ${hapmap}.hmp.txt ${hapmap}