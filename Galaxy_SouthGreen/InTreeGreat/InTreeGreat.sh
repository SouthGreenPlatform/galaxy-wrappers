#!/bin/bash

tool_path=$(dirname $0)

filein=$1
fileout=$2

perl $tool_path/InTreeGreat.pl $filein $fileout

