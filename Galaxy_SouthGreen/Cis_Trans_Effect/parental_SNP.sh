#!/bin/bash


tool_path=$(dirname $0)

filein=$1
fileout=$2
filelog=$3
nameparent1=$4
nameparent2=$5
nb_indv_parent1=$6
nb_indv_parent2=$7
miss_parent1=$8
miss_parent2=$9
num_occ=${10}
min_depth=${11}
min_score=${12}
min_map=${13}

#echo $tool_path/parental_SNP.pl -i $filein -o $fileout -s1 $nameparent1 -s2 $nameparent2 -n1 $nb_indv_parent1 -n2 $nb_indv_parent2 >>$filelog

perl $tool_path/parental_SNP.pl -i $filein -o $fileout -s1 $nameparent1 -s2 $nameparent2 -n1 $nb_indv_parent1 -n2 $nb_indv_parent2 -m1 $miss_parent1 -m2 $miss_parent2 -g $num_occ -p $min_depth -q $min_score -mp $min_map >>$filelog 2>&1 



