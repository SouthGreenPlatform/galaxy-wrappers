#!/bin/bash
argA=$1
argB=$2
shift 2
directory=`dirname $0`
/usr/local/java/jdk6/bin/java -Xmx3000M -cp $directory AlignmentEditor -input ${argA} -output ${argB} -fasta2phylip $*;
