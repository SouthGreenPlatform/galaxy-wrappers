#!/bin/bash
input1=$1
input2=$2
output=$3
name1=$4
name2=$5
logfile=$6

directory=`dirname $0`
 
perl $directory/kallisto2edger.pl -f $input1 -s $input2 -o $output -n $name1 -c $name2 >>$logfile 2>&1

