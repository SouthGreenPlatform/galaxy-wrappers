#!/bin/bash
vcf=$1
gff=$2
type=$3
out=$4
log=$5

directory=`dirname $0`
mkdir tmpdir$$
#cp -rf $input tmpdir$$/input
 
/usr/bin/perl $directory/ConvertVCFinGenomicPositions.pl -v $vcf -g $gff -t $type -o $out >>$log 2>&1



