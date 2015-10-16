#!/bin/bash

# Parameter
input=$1  
keyfile=$2
enzym=$3
max_number_of_good_barcoded=$4
min_tag_number=$5
output=$6 

# Create directory
directory=$(dirname $1) 
output_directory=$6_directory 
input_directory=$6_input
mkdir -p ${output_directory} 
mkdir -p ${input_directory} 
current_directory=$(dirname ${output_directory}) 
base_output_directory=$(basename ${output_directory})
path=$(dirname $0) 
process_name=fq$$
output_name="" 
num_files=0
dos2unix ${keyfile} 2>/dev/null

# For each fastq file in directory, exec FastqToTagCountPlugin
for entry in "$directory"/*gz
do	  
	((num_files++))
	mydir=$(mktemp -p ${input_directory} -d "fastq.XXXXXXXXXX")
	chmod -Rf 777 ${mydir}
	ln -s ${entry} ${mydir}  	
	name=$(basename ${entry})
	#rename .fastq.gz _fastq.gz ${mydir}/*.fastq.gz > /tmp || true
	output_name=${name/_fastq.gz/.cnt}  
	#echo "qsub -V -b y -q web.q -l mem_free=48G -N ${process_name} run_pipeline.pl -Xms512m -Xmx24g -fork1 -FastqToTagCountPlugin  -i ${mydir}/ -o ${output_directory} -k ${keyfile} -e ${enzym} -s ${max_number_of_good_barcoded} -c ${min_tag_number} -endPlugin -runfork1	 "
	qsub -V -b y -q web.q -l mem_free=48G -N ${process_name} run_pipeline.pl -Xms512m -Xmx24g -fork1 -FastqToTagCountPlugin  -i ${mydir}/ -o ${output_directory} -k ${keyfile} -e ${enzym} -s ${max_number_of_good_barcoded} -c ${min_tag_number} -endPlugin -runfork1	 
done  

cnt=$(qstat | grep -c ${process_name})


# Check if all jobs are complete
while [ ${cnt} -gt 0 ]
do  
	cnt=$(qstat | grep -c ${process_name})
  	sleep 10
done 
sleep 10

# If more than one fastq in directory, Merge output
if [ ${num_files} -gt 1 ]
then
	#echo "qsub -V -b y -q web.q -l mem_free=48G -sync yes run_pipeline.pl -Xms512m -Xmx48g -fork1 -MergeMultipleTagCountPlugin -i ${output_directory} -c ${min_tag_number} -o ${output} -endPlugin -runfork1"
	qsub -V -b y -q web.q -l mem_free=48G -sync yes run_pipeline.pl -Xms512m -Xmx24g -fork1 -MergeMultipleTagCountPlugin -i ${output_directory} -c ${min_tag_number} -o ${output} -endPlugin -runfork1
	`cd ${current_directory};rm -Rf ${base_output_directory}` 
else 
	`mv ${output_directory}/${output_name} ${output}`
fi
