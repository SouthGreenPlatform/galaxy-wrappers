#!/bin/bash


vcf_input=$1
missingness=$2
ploidy=$3
vcf_output=$4



mkdir tmpdir$$

cp -rf $vcf_input tmpdir$$/input.vcf


/usr/bin/perl $HOME/galaxy_dist/tools/filters/discard_individuals_based_on_missingness.pl -i tmpdir$$/input.vcf -m $missingness -p $ploidy -o $vcf_output

rm -rf tmpdir$$

