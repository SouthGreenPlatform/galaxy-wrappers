#!/bin/bash

tool_path=$(dirname $0)

filein=$1
filein2=$2
groups=$3
fileout=$4
dotfile=$5
cytoscape_html=$6
logfile=$7


perl $tool_path/Haplophyle.pl --input $filein --groups $groups --stats $filein2 --dot $dotfile --out $fileout --html $cytoscape_html >>$logfile 2>&1



