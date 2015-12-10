#!/bin/bash
vcf_input=$1
reference=$2
depth_of_coverage=$3
min_depth=$4
min_freq=$5
results_zip=$6
heterozygotes=$7


directory=`dirname $0`
 
perl $directory/VCFToFastaAlignments.pl -i $vcf_input -r $reference -d $depth_of_coverage -p tmpdir$$ -n reference -c $min_depth -f $min_freq

mv "tmpdir$$.results.zip" $results_zip;
mv "tmpdir$$.heterozygote_positions.xls" $heterozygotes;


