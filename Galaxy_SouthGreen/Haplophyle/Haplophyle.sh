#!/bin/bash

tool_path=$(dirname $0)

filein=$1
filein2=$2
fileout=$3
dotfile=$4
cytoscape_html=$5
logfile=$6
groups=$7

perl $tool_path/Haplophyle.pl --input $filein --groups $groups --stats $filein2 --dot $dotfile --out $fileout --html $cytoscape_html --tool_path $tool_path >>$logfile 2>&1




