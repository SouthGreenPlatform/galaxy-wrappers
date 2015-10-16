#!/bin/bash

hapmap=$1  
pedf=$2
sC=$3
eC=$4
mxE=$5
mnD=$6
mnPLD=$7
kpUT=$8
output=$9
directory=$(dirname $1)
output_directory=${output}_correction
mkdir -p ${output_directory}

temp_out=" -o ${output_directory}/Merge-hapmap_chr+.hmp.txt"
base_directory=$(basename ${output_directory})
current_directory=$(dirname ${output_directory})
is_kpUT=""
if [ ${kpUT} == true ]
then
	is_kpUT=" -kpUT "
fi

`tar xzf $hapmap`
input_directory=${hapmap}_final/"Final-hapmap_chr+.hmp.txt"
input=" -hmpFile ${input_directory} "

process_name=fq$$
path=$(dirname $0)

for i in $(seq ${sC} ${eC})
do 
	#echo "qsub -b y -q normal.q -N ${process_name} run_pipeline.pl -fork1 -BiParentalErrorCorrectionPlugin ${input}  ${temp_out} -pedF ${pedf} -mxE ${mxE} -mnD ${mnD} -${is_kpUT}  -sC ${i} -eC ${i}  -endPlugin -runfork1 2>&1"
	qsub -b y -q normal.q -N ${process_name} run_pipeline.pl -fork1 -BiParentalErrorCorrectionPlugin ${input} ${temp_out} -pedF ${pedf} -mxE ${mxE} -mnD ${mnD} -${is_kpUT}   -sC ${i} -eC ${i}  -endPlugin -runfork1 2>&1
done
cnt=$(qstat | grep -c ${process_name})
while [ ${cnt} -gt 0 ]
do  
	cnt=$(qstat | grep -c ${process_name})
  	sleep 10
done 
sleep 10
`cd ${current_directory};tar czf ${output} ${base_directory};rm -Rf ${base_directory}`

