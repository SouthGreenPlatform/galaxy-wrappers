#!/bin/bash
scaffolds_fasta=$1
order_tabular=$2
species=$3
number_N=$4
agp=$5
pseudomolecules=$6
gff_pseudo=$7
gff_scaffolds=$8

directory=`dirname $0`

/usr/bin/perl $directory/CreatePseudochromosomes.pl -f $scaffolds_fasta -o $order_tabular -e $species -n $number_N -a $agp -c $pseudomolecules -p $gff_pseudo -s $gff_scaffolds

