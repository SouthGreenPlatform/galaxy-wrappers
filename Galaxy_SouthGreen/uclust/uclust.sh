#!/bin/bash
input=$1
path=$2
output=$3
outputfinal=$4
shift 4
$path/uclust3.0.617_i86linux32 --sort $input --output ${input}_sorted --quiet;
$path/uclust3.0.617_i86linux32 --input ${input}_sorted --uc $output --quiet $*; 
java -cp $path/ ParseClustersGalaxy -sequences $input -clusters $output -output $outputfinal