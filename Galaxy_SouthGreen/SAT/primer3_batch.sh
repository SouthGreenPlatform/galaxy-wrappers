#!/bin/bash

args=("$@")
directory=`dirname $0`

perl $directory/primer3_batch.pl -i $1 -l $2 -o $3 -t $4,$5,$6 -s $7,$8,$9 -g ${args[9]},${args[10]}
