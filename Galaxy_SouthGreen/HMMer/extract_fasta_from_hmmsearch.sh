#!/bin/bash
hmmer_outfile=$1
reference_fastafile=$2
threshold=$3
value=$4
output_file=$5
$HOME/galaxy/tools/SouthGreen/HMMer/extract_fasta_from_hmmsearch.pl -f $reference_fastafile -h $hmmer_outfile -t $threshold -v $value  > $output_file;
