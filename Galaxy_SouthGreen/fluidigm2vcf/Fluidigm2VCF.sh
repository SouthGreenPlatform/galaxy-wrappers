#!/bin/bash
input=$1
output=$2
positions=$3

directory=`dirname $0`
mkdir tmpdir$$

 
perl $directory/Fluidigm2VCF.pl -i $input -o $output -s $positions


