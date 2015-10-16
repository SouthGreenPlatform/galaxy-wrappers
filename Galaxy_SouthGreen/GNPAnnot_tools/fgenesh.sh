#!/bin/bash
input_fasta=$1
models=$2

/usr/local/bioinfo/galaxy/galaxy_dist/tools/prediction/fgenesh.pl $input_fasta $models $3;
mv "$1.raw.fg" $3;
