import sys
import os
import re

def get_field_samples_options(dataset):
	options = []
	line=os.popen("head -20000 %s | grep CHROM"%dataset.file_name).read()[:-1].split('\t')
	index=line.index('FORMAT')
	for opt in line[index+1:] :
		options.append((opt,opt, True))
	return options


