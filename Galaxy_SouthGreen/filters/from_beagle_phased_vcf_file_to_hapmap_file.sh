#!/bin/bash


vcf_input=$1
hapmap_output=$2
del_output=$3

mkdir tmpdir$$

cp -rf $vcf_input tmpdir$$/input.vcf


/usr/bin/perl $HOME/galaxy_dist/tools/filters/from_beagle_phased_vcf_file_to_hapmap_file.pl -i tmpdir$$/input.vcf -o $hapmap_output -d $del_output

rm -rf tmpdir$$

