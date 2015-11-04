#!/bin/bash
input_fam=$1
input_ped=$2
input_simfam=$3
input_simped=$4
fam_suffix=.tfam
ped_suffix=.tped

cp ${input_fam} ${input_fam}${fam_suffix};
cp ${input_ped} ${input_fam}${ped_suffix};
cp ${input_simfam} ${input_simfam}${fam_suffix};
cp ${input_simped} ${input_simfam}${ped_suffix};

directory=`dirname $0`

shift 4

echo "fastlmmc -tfile ${input_fam} -tfileSim ${input_simfam} $* > ${directory}/fastlmm_so.txt &> ${directory}/fastlmm_se.txt;"

fastlmmc -tfile ${input_fam} -tfileSim ${input_simfam} $* > ${directory}"/fastlmm_so.txt" &> ${directory}"/fastlmm_se.txt";

