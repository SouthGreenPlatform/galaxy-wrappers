#!/bin/bash

tool_path=$(dirname $0)

counts=$1
modules=$2
pdf=$3
minsize=$4
threshold=$5
logfile=$6

perl $tool_path/WGCNA.pl -c $counts -m $modules -f $pdf -s $minsize -o $threshold -p $tool_path >>$logfile 2>&1



