#!/bin/bash

# Input parameter
file=$1
type=$2
misMat=$3
callHets=$4 
kpUnmergDups=$5
sC=$6 
eC=$7
maxAlleleVCF=$8
output=$9
snpLog=${10} 
customSNPlog=${11}
is_callHets=""
is_kpUnmergDups=""

# Create directory
dirfile=$(dirname ${file}) 
`cd ${dirfile};tar xzf  $file`
output_directory=${output}_merge
mkdir -p ${output_directory}
file_name=$(basename ${file}) 
base_directory=$(basename ${output_directory})
current_directory=$(dirname ${output_directory})
hapmap_directory=${dirfile}/${file_name}_hapmap
temp_out=""
input=""
 
# Set parameter
if [ ${callHets} == true ]
then
	is_callHets=" -callHets "
fi 

if [ ${kpUnmergDups} == true ]
then
	is_kpUnmergDups=" -kpUnmergDups "
fi 

if [ ${type} == "hapmap" ]
then
	input_directory=${dirfile}/${file_name}_hapmap/"hapmap_chr+.hmp.txt"
	input=" -hmp	${input_directory} "
	temp_out=" -o ${output_directory}/Merge-hapmap_chr+.hmp.txt"
else 
	input_directory=${dirfile}/${file_name}_hapmap/"hapmap_chr+.vcf"
	input=" -vcf ${input_directory} "
	temp_out=" -o ${output_directory}/Merge-hapmap_chr+.vcf"
fi 


# For each chromosome, exec MergeDuplicateSNPsPlugin
process_name=fq$$
for i in $(seq ${sC} ${eC})
do 
	#echo "qsub -V -b y -q web.q -l mem_free=48G -N ${process_name} run_pipeline.pl -Xms512m -Xmx24g  -fork1 -MergeDuplicateSNPsPlugin ${input} ${is_callHets} ${is_kpUnmergDups} ${temp_out} -misMat ${misMat} -snpLog ${snpLog} -sC ${i} -eC ${i}  -endPlugin -runfork1 2>&1"
	qsub -V -b y -q web.q -l mem_free=48G -N ${process_name} run_pipeline.pl -Xms512m -Xmx24g  -fork1 -MergeDuplicateSNPsPlugin ${input} ${is_callHets} ${is_kpUnmergDups} ${temp_out} -misMat ${misMat} -snpLog ${snpLog} -sC ${i} -eC ${i}  -endPlugin -runfork1 2>&1
done


# Check if all jobs are complete
cnt=$(qstat | grep -c ${process_name}) 
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

# Create output
num_chr=0 
for i in $(seq ${sC} ${eC})
do  
	let num_chr=num_chr+1
	if [ $num_chr -eq 1 ]
	then
		#echo "cd ${hapmap_directory};cp hapmap_chr${i}.customSNPLog.txt ${customSNPlog}"
		`cd ${hapmap_directory};cp hapmap_chr${i}.customSNPLog.txt ${customSNPlog}`
	else
		#echo "cd ${hapmap_directory};sed '1d' hapmap_chr${i}.customSNPLog.txt >> ${customSNPlog}"
		`cd ${hapmap_directory};sed '1d' hapmap_chr${i}.customSNPLog.txt >> ${customSNPlog}`
	fi
done

#echo "rm -Rf ${hapmap_directory}"
`rm -Rf ${hapmap_directory}`
#echo "cd ${current_directory};tar czf ${output} ${base_directory};rm -Rf ${base_directory}"
`cd ${current_directory};tar czf ${output} ${base_directory};rm -Rf ${base_directory}` 

