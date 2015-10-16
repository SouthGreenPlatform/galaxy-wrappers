#!/bin/bash

degradome=$1
miRNA=$2
transcripts=$3
output=$4
log=$5

directory=`dirname $0`

export PATH=$PATH:$directory/ViennaRNA-2.1.8/Progs:$directory/CleaveLand4-master/GSTAr_v1-0

perl $directory/CleaveLand4-master/CleaveLand4.pl -e $degradome -u $miRNA -n $transcripts -o outdir -t 1> $output 2> $log 
