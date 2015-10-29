#!/bin/bash

directory=`dirname $0`
#$directory/tgicl/tgicl $1 -q $2 -p $8 -l $9 -O '-p '$8' -o '$9 >>${10} 2>&1
tgicl $1 -q $2 -p $8 -l $9 -O '-p '$8' -o '$9 >>${10} 2>&1

perl $directory/contigtoXML.pl $PWD

echo $PWD >>${10}

cat asm_*/ACE >$3
cat asm_*/contigs >$4

#/usr/local/bioinfo/galaxy/galaxy_dist/tools/ESTtik/tgicl/bin/formatdb -i $1 -p F
#/usr/local/bioinfo/galaxy/galaxy_dist/tools/ESTtik/fastacmd -d $1 -i $1.singletons >$5

cp -rf asm_1/singlets $5 2>tgicl.err
cat *.singletons >>$5 2>>tgicl.err
cat $4 >$6
cat $5 >>$6
cp -rf asm_1/align $7



