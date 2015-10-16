#!/bin/bash

input=$1 
output=$2 
dirfile=$(dirname ${input}) 
`cd ${dirfile};tar xzf  $input`
file_name=$(basename ${input})
input_directory=${dirfile}/${file_name}_merge
 
vcf-concat ${input_directory}/*vcf | bgzip > ${output}
rm -Rf ${input_directory}
