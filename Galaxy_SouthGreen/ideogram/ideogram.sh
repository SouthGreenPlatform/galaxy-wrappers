#!/bin/bash

tool_path=$(dirname $0)
chromosomelength=$1
annotation=$2
chrfile=${chromosomelength##*/}
annotile=${annotation##*/}
ploidy=$3
fileout=$4

perl $tool_path/ideogram.pl $chromosomelength $annotation $ploidy $fileout

