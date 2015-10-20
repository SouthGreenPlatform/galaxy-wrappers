#!/bin/bash


vcf_input=$1
reference=$2
vcf_output=$3


mkdir tmpdir$$

cp -rf $vcf_input tmpdir$$/input.vcf
cp -rf $reference tmpdir$$/reference.fasta


/usr/bin/perl $HOME/galaxy_dist/tools/filters/cornell_vcf_to_standard_vcf.pl -i tmpdir$$/input.vcf -r tmpdir$$/reference.fasta -o $vcf_output

rm -rf tmpdir$$

