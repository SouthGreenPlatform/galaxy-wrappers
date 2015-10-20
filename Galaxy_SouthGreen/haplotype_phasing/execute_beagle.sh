#!/bin/bash


vcf_input=$1
phase=$2
impute=$3
out_prefix="out"



mkdir tmpdir$$

cp -rf $vcf_input tmpdir$$/input.vcf


/usr/local/jdk1.7.0_40/bin/java -Xmx2g -jar /usr/local/bioinfo/galaxy_dev/galaxy_dist/tools/GWAS_analysis/beagle/beagle.r1399.jar gt=$vcf_input phase-its=$phase impute-its=$impute out=$out_prefix

gunzip $out_prefix.vcf.gz

rm -rf tmpdir$$

