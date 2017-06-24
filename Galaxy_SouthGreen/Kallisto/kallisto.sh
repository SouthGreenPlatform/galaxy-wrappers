#!/bin/bash
reference=$1
forward=$2
reverse=$3
abundance=$4
logfile=$5

directory=`dirname $0`
 
$directory/kallisto_linux-v0.43.1/kallisto index -i $reference.index $reference >>$logfile 2>&1
$directory/kallisto_linux-v0.43.1/kallisto quant -i $reference.index -o output$$ $forward $reverse >>$logfile 2>&1 

mv output$$/abundance.tsv $abundance;
