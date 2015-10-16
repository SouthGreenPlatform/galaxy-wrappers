#!/bin/bash
input=$1
tool_directory=$2
model=$3
output_mst=$4
output_mad=$5
 

if [ "$model" = "rice" ]
then
donors=${tool_directory}/model.rice.donors.w;
acceptors=${tool_directory}/model.rice.acceptors.w;
starts=${tool_directory}/model.rice.starts.w;
options="-p  -P1_gt 80 80  -P3_gt 40 100  -C4_gt 80 60  -C3_gt 5 5  -P2_gt 40 120 -C6_gt 80 100  -RF_gt 5 5 -P1_ag 80 80  -C5_ag 80 80  -P3_ag 60 80  -RF_ag 80 60  -P1_ag 80 80 -C3_ag 60 60  -P2_ag 80 100  -C6_ag 100 80 -P2_atg 120 120  -P1_atg 100 80  -C6_atg 120 100 ";
else
donors=${tool_directory}/model.medicago.donors.w;
acceptors=${tool_directory}/model.medicago.acceptors.w;
starts=${tool_directory}/model.medicago.starts.w;
options="-p -P1_ag 70 65 -C4_ag 50 75 -C6_ag 75 80 -RF_ag 50 80 -P3_ag 70 55 -P2_ag 80 70 -sig_ag -2.388988 -0.065396 -P1_gt 80 85 -RF_gt 65 65 -P3_gt 75 65 -C5_gt 90 65 -sig_gt -2.320267 -0.068317 -RF_atg 65 75 -P2_atg 65 110 -C3_atg 115 90 -C4_atg 115 85 -P1_atg 50 110 -P3_atg 55 120 -sig_atg -2.616002 -0.213646";
fi



${tool_directory}/splicemachine $options -model_gt $donors -model_ag $acceptors -model_atg $starts $1;

mv "$1.spliceMSt" $4;
mv "$1.spliceMAD" $5;


