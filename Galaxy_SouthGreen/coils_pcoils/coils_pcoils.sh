#!/bin/bash
fasta=$1
window_length=$2
output=$3
log=$4
threshold=$5



directory=`dirname $0`

COILSDIR=$directory/PCOILS_v1.0.1
export COILSDIR

/usr/bin/perl $directory/coils_pcoils_batch.pl -f $fasta -w $window_length -t $threshold -o $output -p $directory >>$log 2>&1



