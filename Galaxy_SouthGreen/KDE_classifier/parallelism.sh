#!/bin/bash
qsub -N FM36 -V -q bigmem.q -pe parallel_smp 10 -l mem_free=20G -b y /usr/local/bioinfo/python/2.7.9_build2/bin/python $1 
--CHR $2
--geno $3
--fam $4
--bim $5
--ref $6
--admx $7
--proc 10;

