#!/bin/bash
ace_input=$1
min_depth=$2
hetero_depth=$3
occ_min=$4
freq_min=$5
alignments_zip=$6
heterozygotes_sites=$7
SNP=$8

directory=`dirname $0`

perl $directory/AceToFastaAlignments.pl -i $ace_input -m $min_depth -h $hetero_depth -o $occ_min -f $freq_min -p tmpdir$$

mv "tmpdir$$.set0.fasta_alignments.zip" $alignments_zip
mv "tmpdir$$.heterozygous_sites" $heterozygotes_sites
mv "tmpdir$$.SNP" $SNP

