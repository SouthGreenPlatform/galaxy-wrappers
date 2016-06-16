import sys
import os
import re

def get_field_samples_options(dataset):
	options = []
	line=os.popen("grep '#CHROM' %s"%dataset.file_name).read()[:-1].split('\t')
	index=line.index('FORMAT')
	for opt in line[index+1:] :
		options.append((opt,opt, True))
	return options

def get_field_chrs_options(dataset):
        options = []
        chrs=os.popen("grep '##contig' %s"%dataset.file_name).read()[:-1].split('\n')
        for line in chrs:
		opt=re.search('^##contig=<ID=(\w+),length=',line).group(1)
                options.append((opt,opt, True))
        return options

