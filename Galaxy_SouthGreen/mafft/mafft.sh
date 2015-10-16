#!/bin/bash


directory=`dirname $0`
output_tree=$1
input=$2
tree_suffixe=.tree
tmp_suffixe=.tmp
shift 2
mafft --treeout $*;
sed 's/$/;/g' ${input}${tree_suffixe} > ${output_tree}
rm ${input}${tree_suffixe}


