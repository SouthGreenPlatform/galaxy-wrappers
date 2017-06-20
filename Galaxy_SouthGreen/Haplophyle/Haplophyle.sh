#!/bin/bash

tool_path=$(dirname $0)

filein=$1
fileout=$2
dotfile=$3
cytoscape_html=$4
logfile=$5


perl $tool_path/Haplophyle.pl --input $filein --dot $dotfile --out $fileout --html $cytoscape_html >>$logfile 2>&1



