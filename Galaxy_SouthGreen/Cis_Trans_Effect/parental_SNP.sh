#!/bin/bash


tool_path=$(dirname $0)

filein=$1
fileout=$2
filelog=$3
nameparent1=$4
nameparent2=$5
nb_indv_parent1=$6
nb_indv_parent2=$7


#echo $tool_path/parental_SNP.pl -i $filein -o $fileout -s1 $nameparent1 -s2 $nameparent2 -n1 $nb_indv_parent1 -n2 $nb_indv_parent2 >>$filelog

perl $tool_path/parental_SNP.pl -i $filein -o $fileout -s1 $nameparent1 -s2 $nameparent2 -n1 $nb_indv_parent1 -n2 $nb_indv_parent2 >>$filelog 2>&1 



