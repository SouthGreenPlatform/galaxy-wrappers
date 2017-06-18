#!/bin/bash

tool_path=$(dirname $0)

filein=$1
fileout=$2
logfile=$3


perl $tool_path/Haplophyle.pl --input $filein --out $fileout >>$logfile 2>&1



