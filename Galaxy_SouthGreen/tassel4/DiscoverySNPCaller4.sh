#!/bin/bash

# Input Parameter
input=$1
useTBTByte=$2
topm=$3
mUpd=$4
vcf=$5
mxSites=$6
mnF=$7
mnMAF=$8 
mnMAC=$9
mnLCov=${10}
errRate=${11}
sC=${12}
eC=${13}
inclRare=${14}
inclGaps=${15} 
callBiSNPsWGap=${16}
database=${17}
output_hapmap=${18} 
 
# Create directory
directory=$(dirname $1)
hapmap_directory=${output_hapmap}_hapmap
current_directory=$(dirname ${output_hapmap})
mkdir -p ${hapmap_directory}
base_hapmap_directory=$(basename ${hapmap_directory})
hapmap_file=${hapmap_directory}/"hapmap_chr+.hmp.txt" 


# Eval several parameter
is_tbtbyte=""
is_vcf=""
is_inclRare=""
is_inclGaps=""
is_callbisnp=""
if [ $useTBTByte == true ]
then
	is_tbtbyte=" -y "
fi 

if [ $vcf == true ]
then
	is_vcf=" -vcf " 
fi
if [ $inclRare == true ]
then
	is_inclRare=" -inclRare "
fi
if [ $inclGaps == true ]
then
	is_inclGaps=" -inclGaps "
fi
if [ $callBiSNPsWGap == true ]
then
	is_callbisnp=" -callBiSNPsWGap "
fi
 
process_name=fa$$ 

# For each chromosome, exec DiscoverySNPCaller
for i in $(seq ${sC} ${eC})
do  
	#echo "qsub -b y -q web.q -l mem_free=48G -N ${process_name} -V run_pipeline.pl -Xms512m -Xmx36g -fork1 -DiscoverySNPCallerPlugin -i ${input} ${is_callbisnp} ${is_inclRare} ${is_inclGaps} ${is_vcf} ${is_tbtbyte} -m ${topm} -o ${hapmap_file} -ref ${database} -sC ${i} -eC ${i} -mnLCov ${mnLCov} -mnMAC ${mnMAC} -errRate ${errRate} -mnMAF ${mnMAF}  -endPlugin -runfork1 "
	qsub -b y -q web.q -l mem_free=48G -N ${process_name} -V run_pipeline.pl -Xms512m -Xmx36g -fork1 -DiscoverySNPCallerPlugin -i ${input} ${is_callbisnp} ${is_inclRare} ${is_inclGaps} ${is_vcf} ${is_tbtbyte} -m ${topm} -o ${hapmap_file} -ref ${database} -sC ${i} -eC ${i} -mnLCov ${mnLCov} -mnMAC ${mnMAC} -errRate ${errRate} -mnMAF ${mnMAF}  -endPlugin -runfork1 

done

# Check if all jobs are complete
cnt=$(qstat | grep -c ${process_name})

while [ ${cnt} -gt 0 ]
do  
	cnt=$(qstat | grep -c ${process_name})
  	sleep 10
done 
sleep 10
if [ $vcf == true ]
echo $vcf
then
	#echo "cd ${hapmap_directory};sed -i 's/####/##/' *.vcf"
	cd ${hapmap_directory};sed -i 's/####/##/' *.vcf
fi

sleep 10
#echo "cd ${current_directory};tar czf ${output_hapmap} ${base_hapmap_directory}" 
`cd ${current_directory};tar czf ${output_hapmap} ${base_hapmap_directory}; rm -Rf ${base_hapmap_directory}` 
