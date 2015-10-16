#!/bin/bash
input=$1
species=$2
bank=$3
output=$4
if [[ $bank == "/bank/est/MS_mrnas.fsa" ]]
	then suffix="MS_mrnas"
elif  [[ $bank == "/bank/est/OS_embl_mrnas.fsa" ]]
	then suffix="OS_mrnas"
elif [[ $bank == "/bank/est/WH_embl_mrna.fsa" ]]
	then suffix="WH_embl"
elif [[ $bank == "/bank/est/ZM_embl_mrnas.fsa" ]]
	then suffix="ZM_mrnas"
elif [[ $bank == "/bank/est/SB_embl_mrnas.fsa" ]]
	then suffix="SB_mrnas"
elif [[ $bank == "/bank/est/SC_embl_mrnas.fsa" ]]
	then suffix="SC_mrnas"
elif [[ $bank == "/bank/est/RU_mrnas.fsa" ]]
	then suffix="RU_mrnas"
elif [[ $bank == "/home/sidibebocs/work/Coffea/454/clean/Ca_454.fna" ]]
	then suffix="CC_mrnas"
elif [[ $bank == "/bank/est/TE_mrnas.fsa" ]]
	then suffix="TE_mrnas"
else suffix="Other"
fi

if [ -e "$input.tis" ]
  then
 	gth -genomic $input -cdna $bank  -o $output -bssm $species -force yes -noautoindex
  else  
	gth -genomic $input -cdna $bank  -o $output -bssm $species -force yes -createindicesonly	
	gth -genomic $input -cdna $bank  -o $output -bssm $species -force yes
fi
