#!/bin/bash
input=$1
output_aln=$2
output_htm=$3
ba=$4
bb=$5
aln_suffix=-gb
htm_suffix=-gb.htm
shift 5

directory=`dirname $0`
/usr/local/java/jdk6/bin/java -mx4000m -cp $directory/ GetNbSequences $input | \
while read nbseq
do
	$directory/Gblocks $input -b1=$((nbseq * ba / 100 + 1)) -b2=$((nbseq * bb / 100 + 1)) $*;
	mv ${input}${aln_suffix} ${output_aln};
	mv ${input}${htm_suffix} ${output_htm};	
done

