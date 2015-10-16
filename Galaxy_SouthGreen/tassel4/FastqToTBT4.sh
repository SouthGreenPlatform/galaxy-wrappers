#!/bin/bash
 
# Input parameter
filename=$1  
keyfile=$2
enzym=$3 
topm=$4 
min_taxa_number=$5
max_number_of_tags=$6
merge_tag_count=$7
output=$8
logfile=$9 

# Create directory
directory=$(dirname $1)
output_directory=${output}_directory
input_directory=${output}_input
mkdir -p ${output_directory} 
mkdir -p ${input_directory} 

dos2unix ${keyfile} 2>/dev/null
 
is_x="" 
if [ ${merge_tag_count} == true ]
then
	is_x=" -x "
fi
 

# For each fastq file in directory, exec FastqToTBTPlugin
process_name=fq$$
num_files=0
for entry in "$directory"/*gz
do	  
	((num_files++))
	mydir=$(mktemp -p ${input_directory} -d "fastq.XXXXXXXXXX")
	chmod -Rf 777 ${mydir}
	ln -s ${entry} ${mydir}  
	#echo "qsub -V -b y -q web.q -l mem_free=32G -N ${process_name} run_pipeline.pl -Xms512m -Xmx24g -fork1 -FastqToTBTPlugin  -i ${mydir} -o ${output_directory} -k ${keyfile} -e ${enzym} -y -m ${topm} -c ${min_taxa_number} -endPlugin -runfork1  > ${logfile} 2>&1" 
	qsub -V -b y -q web.q -l mem_free=32G -N ${process_name} run_pipeline.pl -Xms512m -Xmx24g -fork1 -FastqToTBTPlugin  -i ${mydir} -o ${output_directory} -k ${keyfile} -e ${enzym} -y -m ${topm} -c ${min_taxa_number} -endPlugin -runfork1  2>&1 
done  


# Check if all jobs are complete
cnt=$(qstat | grep -c ${process_name})
echo ${cnt}

while [ ${cnt} -gt 0 ]
do  
	cnt=$(qstat | grep -c ${process_name})
  	sleep 10
done 
sleep 10
 
 
# If more than one fastq in directory, Merge output
temp_file=$(mktemp -p ${input_directory} "merge.XXXXXXXXXX.tbt.byte")
#echo "qsub -V -b y -q web.q -l mem_free=32G -N ${process_name} -sync yes run_pipeline.pl -Xms512m -Xmx24g  -fork1 -MergeTagsByTaxaFilesPlugin -i ${output_directory} -s ${max_number_of_tags} ${is_x} -o ${temp_file} -endPlugin -runfork1"
qsub -V -b y -q web.q -l mem_free=32G -N ${process_name} -sync yes run_pipeline.pl -Xms512m -Xmx24g  -fork1 -MergeTagsByTaxaFilesPlugin -i ${output_directory} -s ${max_number_of_tags} ${is_x} -o ${temp_file} -endPlugin -runfork1 >> ${logfile} 2>&1
#echo "mv ${temp_file} ${output}"
mv ${temp_file} ${output}
