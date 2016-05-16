#!/bin/sh

VERT="\\033[1;32m" NORMAL="\\033[0;39m" 
ROUGE="\\033[1;31m" JAUNE="\\033[1;33m"

cd code

make clean
make

cd ..

if [ ! -e "code/lib/lapack.a" ]; then
  	echo "$ROUGE" "\n ERROR: an error occured during the compilation of lapack library." "$NORMAL"
        exit 1
elif [ ! -e "bin/sNMF" ]; then
	echo "$ROUGE" "\n ERROR: an error occured during the compilation of sNMF command-line program." "$NORMAL"
	exit 1
elif [ ! -e "bin/crossEntropy" ]; then
  	echo "$ROUGE" "\n ERROR: an error occured during the compilation of crossEntropy command-line program." "$NORMAL"
        exit 1
elif [ ! -e "bin/createDataSet" ]; then
  	echo "$ROUGE" "\n ERROR: an error occured during the compilation of createDataSet command-line program." "$NORMAL"
        exit 1
elif [ ! -e "bin/vcf2geno" ]; then
  	echo "$ROUGE" "\n ERROR: an error occured during the compilation of vcf2geno command-line program." "$NORMAL"
        exit 1
elif [ ! -e "bin/ancestrymap2geno" ]; then
  	echo "$ROUGE" "\n ERROR: an error occured during the compilation of ancestrymap2geno command-line program." "$NORMAL"
        exit 1
elif [ ! -e "bin/ped2geno" ]; then
  	echo "$ROUGE" "\n ERROR: an error occured during the compilation of ped2geno command-line program." "$NORMAL"
        exit 1
elif [ ! -e "bin/lfmm2geno" ]; then
  	echo "$ROUGE" "\n ERROR: an error occured during the compilation of lfmm2geno command-line program." "$NORMAL"
        exit 1
elif [ ! -e "bin/geno2lfmm" ]; then
  	echo "$ROUGE" "\n ERROR: an error occured during the compilation of geno2lfmm command-line program." "$NORMAL"
        exit 1
fi

echo "$VERT" "\n SUCCESS: sNMF command-line program was compiled without error." "$NORMAL"
exit 0

