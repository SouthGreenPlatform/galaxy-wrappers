#!/usr/bin/env python
# -*- coding: utf8 -*-


from multiprocessing import Process, Manager
import os, sys, random, re
from Bio import SeqIO

input_file = sys.argv[1]
list1=[]
spec_codes=[]

fasta_sequences = SeqIO.parse(open(input_file),'fasta')
for fasta in fasta_sequences:
    list1.append(fasta.id)


output_fna=sys.argv[2]
output_bed=sys.argv[3]
flank=sys.argv[4]
begin=sys.argv[5]
end=sys.argv[6]

for i in list1:
	toto = i.split('_')
	chaine = toto[len(toto)-1]
	expression = r"[A-Z]{5}"
	if re.match(expression, chaine) is not None:
		spec_codes.append(toto[len(toto)-1])
	else:
		print ("wrong id format")

myset = set(spec_codes)
print(myset)

	
#path_file=open('/usr/local/bioinfo/galaxy_dev/galaxy_dist/tool-data/gff2fasta.loc','r')
path_file=open('/home/galaxydev/galaxy/tool-data/SouthGreen/gff2fasta.loc','r')
spec=path_file.readlines()
path_file.close()
bases_gff=[]
bases_fna=[]
for line in spec :
	infos=line.split('\t')
	if infos[0] in spec_codes:
		bases_gff.append(infos[2].rstrip())
		bases_fna.append(infos[3].rstrip())
		print(infos[2].rstrip()+"\t"+infos[3].rstrip())

curdir=os.getcwd()
    	
numb=str(random.randint(1,100000))
tmpfoldname=curdir+"/"+numb
os.system("mkdir "+tmpfoldname)
os.system("mkdir "+tmpfoldname+"/logs")
os.system("mkdir "+tmpfoldname+"/outputs")

#os.chdir(tmpfoldname+"/outputs")
#print os.getcwd()
#os.chdir("../../")

# dirlist = os.listdir(tmpfoldname+"/outputs")
# 
# from pprint import pprint
# pprint(dirlist)

i=1
for path in bases_gff:
	os.system("ln -s "+ path +" "+tmpfoldname+"/gff"+str(i))
	#print path
	i=i+1
	
j=1	
for path in bases_fna:
	os.system("ln -s "+ path +" "+tmpfoldname+"/fna"+str(j))
	#print path
	j=j+1


#marmadais:
#jobarray="#!/bin/bash\n\n#$ -N gff2fnaclust\n#$ -wd "+tmpfoldname+"/\n#$ -e "+tmpfoldname+"/logs/\n#$ -o "+tmpfoldname+"/logs/\n#$ -q bioinfo.q\n#$ -t 1-"+str(i-1)+"\n#$ -tc "+str(i-1)+"\n#$ -S /bin/bash\n#$ -b y\n#$ -V\n\nperl /usr/local/bioinfo/galaxy_dev/galaxy_dist/tools/gff2fna/multi_gff2fna.pl "+tmpfoldname+"/gff${SGE_TASK_ID} "+tmpfoldname+"/fna${SGE_TASK_ID} "+tmpfoldname+"/outputs/output_fna${SGE_TASK_ID} "+tmpfoldname+"/outputs/output_bed${SGE_TASK_ID} "+input_file+" "+flank+" -begin="+begin+" -end="+end+"\n"

#print (jobarray)

#cc2-login:
jobarray="#!/bin/bash\n\n#$ -N gff2fnaclust\n#$ -wd "+tmpfoldname+"/\n#$ -e "+tmpfoldname+"/logs/\n#$ -o "+tmpfoldname+"/logs/\n#$ -q normal.q\n#$ -t 1-"+str(i-1)+"\n#$ -tc "+str(i-1)+"\n#$ -S /bin/bash\n#$ -b y\n#$ -V\n#$ -l h_vmem=8G\n\nperl /home/galaxydev/galaxy/tools/SouthGreen/gff2fna/multi_gff2fna.pl "+tmpfoldname+"/gff${SGE_TASK_ID} "+tmpfoldname+"/fna${SGE_TASK_ID} "+tmpfoldname+"/outputs/output_fna${SGE_TASK_ID} "+tmpfoldname+"/outputs/output_bed${SGE_TASK_ID} "+input_file+" "+flank+" -begin="+begin+" -end="+end+"\n"

array_file=open(tmpfoldname+'/jobs.sge','w')
array_file.write(jobarray)
array_file.close()


# dirlist = os.listdir(tmpfoldname)
# 
# from pprint import pprint
# pprint(dirlist)

os.system("chmod -R 777 "+tmpfoldname)

os.system("qsub -sync y "+tmpfoldname+"/jobs.sge")
#os.chdir(tmpfoldname)


k=1
fo = open(output_fna, "w")
for path in bases_fna:
	file_fna=tmpfoldname+"/outputs/output_fna"+str(k)
	with open(file_fna, 'r') as readfile:
		for line in readfile:
			fo.seek(0, 2)
			fo.writelines( line )
	k=k+1

l=1
fo1 = open(output_bed, "w")
for path in bases_gff:
	file_bed=tmpfoldname+"/outputs/output_bed"+str(l)
	with open(file_bed, 'r') as readfile:
		for line in readfile:
			fo1.seek(0, 2)
			fo1.writelines( line )
	l=l+1

fo.close()
fo1.close()

#os.chdir("../")
os.system("rm -rf "+tmpfoldname)
