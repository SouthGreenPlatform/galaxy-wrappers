#!/bin/bash
directory=`dirname $0`
module load system/java/jdk8
java -cp $directory/ -mx4096m ProteinStructureDraw $*;
