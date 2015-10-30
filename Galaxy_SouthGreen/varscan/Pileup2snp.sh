#!/bin/bash
pileup_input=$1
snp=$2
stats=$3
min_coverage=$4
min_reads2=$5
min_avg_qual=$6
min_var_freq=$7
min_freq_for_hom=$8

directory=`dirname $0`
	
java -Xmx4g -jar $directory/VarScan.v2.2.3.jar pileup2snp $1 --min-coverage $4 --min-reads2 $5 --min-avg-qual $6 --min-var-freq $7 --min-freq-for-hom $8 1>$2 2>$3 




