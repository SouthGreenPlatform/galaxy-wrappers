#!/usr/bin/env python
# -*- coding: utf8 -*-

import sys

input=sys.argv[1]
output= sys.argv[2]
file_to_clean = open(input,'r')
res_file = open(output,'w')

for line in file_to_clean:
    line = line.replace(";",".")
    line = line.replace("/","_")
    array_line=line.split('\t')
    if len(array_line)==8 and array_line[6]!="\"\"":
        print array_line[6]
        res_file.write(line)

file_to_clean.close()
res_file.close()