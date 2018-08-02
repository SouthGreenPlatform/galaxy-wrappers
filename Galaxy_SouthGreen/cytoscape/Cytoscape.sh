#!/bin/bash

tool_path=$(dirname $0)

filein=$1
cytoscape_html=$2
logfile=$3

perl $tool_path/Cytoscape.pl --input $filein --html $cytoscape_html >>$logfile 2>&1




