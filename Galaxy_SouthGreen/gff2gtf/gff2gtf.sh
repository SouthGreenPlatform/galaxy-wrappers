#!/bin/bash
gff3=$1
gtf=$2
log=$3

directory=`dirname $0`
/usr/bin/perl $directory/gff2gtf.pl $gff3 >>$log 2>&1

mv "$1.gtf" $2;

