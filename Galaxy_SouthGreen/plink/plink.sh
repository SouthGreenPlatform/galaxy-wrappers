#!/bin/bash


tool_path=$(dirname $0)

filein=$1
fileout_label=$(date "+%Y%m%d%H%M%S")
fileout=$2
filelog=$3
frequency=$4
max_freq=$5
allow_missing=$6
type=${7}
bound_start=${8}
bound_end=${9}

cp -rf $filein input$$.vcf

if [ "${10}" != "None" ]
then samples="--samples ${10}"
fi

if [ "${11}" != "None" ]
then chromosomes="--chromosomes ${11}"
fi

if [ "$bound_start" -gt "$bound_end" ]
then tmp=$bound_start ; bound_start=$bound_end ; bound_end=$tmp ; echo "Warning : Lower bound must be lower than greater bound!" >&2
fi


export="VCF"

perl $tool_path/Plink.pl --input input$$.vcf --out $fileout_label --export $export --frequency $frequency --max_freq $max_freq --allow_missing $allow_missing --type $type --bounds $bound_start','$bound_end $samples $chromosomes


#echo ${16} >>$fileout_label.log
#echo ${15} >>$fileout_label.log
#echo ${17} >>$fileout_label.log
#echo ${18} >>$fileout_label.log

if [ "$export" = "VCF" ]
then cp  $fileout_label.vcf $fileout ; rm $fileout_label.vcf
else cp  $fileout_label.bed $fileout; cp $fileout_label.bed ${15} ; cp $fileout_label.bim ${18} ;rm $fileout_label.bed $fileout_label.fam $fileout_label.bim
fi

cp $fileout_label.log $filelog
rm $fileout_label.log

