#!/bin/bash

sleep 120 
./structure -m  -e ./structure_analysis.sh -i infile -K 2 -o /work/galaxydev/tmpdir34341/resultat_2_${SGE_TASK_ID}
SEP 2 
./structure -m  -e ./structure_analysis.sh -i infile -K 3 -o /work/galaxydev/tmpdir34341/resultat_3_${SGE_TASK_ID}
SEP 3 
./structure -m  -e ./structure_analysis.sh -i infile -K 4 -o /work/galaxydev/tmpdir34341/resultat_4_${SGE_TASK_ID}
SEP 4 
