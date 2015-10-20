#!/bin/bash
directory=`dirname $0`
java -cp $directory/ -mx4096m ProteinStructureDraw $*;
