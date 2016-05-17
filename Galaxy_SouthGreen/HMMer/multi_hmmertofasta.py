#!/usr/bin/env python
# -*- coding: utf8 -*-


from multiprocessing import Process, Manager
import os, sys, random


hmmer_out = sys.argv[1]
threshold=sys.argv[2]
value=sys.argv[3]
output_file=sys.argv[4]
bases=sys.argv[5:]
    	
curdir=os.getcwd()
numb=str(random.randint(1,100000))
tmpfoldname=curdir+"/"+numb
os.system("mkdir "+tmpfoldname)
os.system("mkdir "+tmpfoldname+"/logs")
os.system("mkdir "+tmpfoldname+"/outputs")


# i=1
# for path in bases:
# 	os.system("ln -s "+ path +" "+tmpfoldname+"/base"+str(i))
# 	print path
# 	i=i+1

os.system("tar -xf "+hmmer_out+" -C "+tmpfoldname+"/")

base_file=open(tmpfoldname+'/bases.txt','r')

i=1
path=base_file.readline()
while(path!=""):
    path=path.replace('\n', '')
    os.system("ln -s "+ path +" "+tmpfoldname+"/base"+str(i))
    print path
#     print "ln -s "+ path +" "+tmpfoldname+"/base"+str(i)
    path=base_file.readline()
    i=i+1

base_file.close()


jobarray="#!/bin/bash\n#$ -N hmmerclust\n#$ -wd "+tmpfoldname+"/\n#$ -e "+tmpfoldname+"/logs/ \n#$ -o "+tmpfoldname+"/outputs/\n#$ -q normal.q\n#$ -t 1-"+str(i-1)+"\n#$ -tc 4\n#$ -b y\n#$ -V\n#$ -S /bin/bash \n/homedir/galaxydev/galaxy/tools/SouthGreen/HMMer/extract_fasta_from_hmmsearch.pl -f "+tmpfoldname+"/base${SGE_TASK_ID} -h "+tmpfoldname+"/outputs/output${SGE_TASK_ID} -t "+threshold+" -v "+value+"  > "+tmpfoldname+"/outs${SGE_TASK_ID}"
array_file=open(tmpfoldname+'/jobs.sge','w')
array_file.write(jobarray)
array_file.close()

# os.system("cp "+tmpfoldname+"/jobs.sge /usr/local/bioinfo/galaxy/galaxy_dist/tools/genfam/testgalaxy.txt")

os.system("chmod +x "+tmpfoldname+"/jobs.sge")
os.system("qsub -sync y "+tmpfoldname+"/jobs.sge")
os.system("cat "+tmpfoldname+"/outs* >> "+output_file) 

#os.system("rm -rf "+tmpfoldname)

