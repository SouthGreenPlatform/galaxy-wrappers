#!/bin/bash
file=$1
type=$2
mnTCov=$3
mnSCov=$4
mnF=$5
mnMAF=$6
mxMAF=$7
hLD=$8
mnR2=$9
mnBonP=${10} 
sC=${11}
eC=${12}
maxAlleleVCF=${13}
output=${14}
snpLog=${15}

dirfile=$(dirname ${file})
#echo "tar xzf $file"
`tar xzf $file`
output_directory=${output}_final
mkdir -p ${output_directory}
base_directory=$(basename ${output_directory})
current_directory=$(dirname ${output_directory})
temp_out=""
input=""
 
directory=$(dirname $1)
if [ ${type} == "hapmap" ]
then
	input_directory=${dirfile}/${file}_merge/"Merge-hapmap_chr+.hmp.txt"
	input=" -hmp ${input_directory} "
	temp_out=" -o ${output_directory}/Final-hapmap_chr+.hmp.txt"
else 
	input_directory=${dirfile}/${file}_merge/"Merge-hapmap_chr+.vcf"
	input=" -vcf ${input_directory} "
	temp_out=" -o ${output_directory}/Merge-hapmap_chr+.hmp.vcf"
fi

path=$(dirname $0)

 
process_name=fq$$
path=$(dirname $0)
for i in $(seq ${sC} ${eC})
do 
	#echo "qsub -b y -q normal.q -N ${process_name} run_pipeline.pl -fork1 -GBSHapMapFiltersPlugin ${input} ${temp_out}  -mnTCov ${mnTCov} -mnSCov ${mnSCov} -mnF ${mnF} -mnMAF ${mnMAF} -mxMAF ${mxMAF} -sC ${i} -eC ${i} -snpLog ${snpLog} -endPlugin -runfork1 2>&1"
	qsub -b y -q normal.q -N ${process_name} run_pipeline.pl -fork1 -GBSHapMapFiltersPlugin ${input} ${temp_out}  -mnTCov ${mnTCov} -mnSCov ${mnSCov} -mnF ${mnF} -mnMAF ${mnMAF} -mxMAF ${mxMAF} -sC ${i} -eC ${i} -snpLog ${snpLog} -endPlugin -runfork1 2>&1
done

cnt=$(qstat | grep -c ${process_name})
echo ${cnt}

while [ ${cnt} -gt 0 ]
do  
	cnt=$(qstat | grep -c ${process_name})
  	sleep 10
done 
sleep 10

if [ ${type} == "vcf" ]
then
	#echo "cd ${output_directory};sed -i 's/####/##/' *.vcf"
	`cd ${output_directory};sed -i 's/####/##/' *.vcf`
fi

sleep 10
#echo "cd ${current_directory};tar czf ${output} ${base_directory};rm -Rf ${base_directory}"
`cd ${current_directory};tar czf ${output} ${base_directory};rm -Rf ${base_directory}` 