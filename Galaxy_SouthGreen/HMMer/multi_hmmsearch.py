#!/usr/bin/env python
# -*- coding: utf8 -*-

from multiprocessing import Process, Manager
import os, sys, random, datetime


profil_path = sys.argv[1]
output_path=sys.argv[2]
if sys.argv[3]=='base_list':
    bases=sys.argv[4:]
elif sys.argv[3]=='file':
    filename=sys.argv[4]
    list_file=open(filename,'r')
    fullchain=list_file.read()
    list_file.close()
    spec_codes=fullchain.rstrip().split(',')
    print (spec_codes[0]+"\n"+spec_codes[1]+"\n"+spec_codes[2]+"\n"+spec_codes[3]+"\n")
    path_file=open('/homedir/galaxydev/galaxy/tool-data/HMMdb_p_annot.loc','r')
    spec=path_file.readlines()
    print (spec[0])
    path_file.close()
    bases=[]
    for line in spec :
        infos=line.split('\t')
        if infos[0] in spec_codes:
            bases.append(infos[2].rstrip())
            print(infos[2].rstrip())
    
    	
curdir=os.getcwd()
random.seed()
numb=str(random.randint(1,1000000))
tmpfoldname=curdir+"/"+numb
print tmpfoldname
os.system("mkdir "+tmpfoldname)
os.system("mkdir "+tmpfoldname+"/logs")
os.system("mkdir "+tmpfoldname+"/outputs")

base_file=open(tmpfoldname+'/bases.txt','w')
log=open(tmpfoldname+'/log.txt','w')

i=1
for path in bases:
	os.system("ln -s "+ path +" "+tmpfoldname+"/base"+str(i))
	print path
	base_file.write(path+"\n")
	i=i+1

base_file.close()


jobarray="#!/bin/bash\n#$ -N hmmerclust\n#$ -wd "+tmpfoldname+"/\n#$ -e "+tmpfoldname+"/logs/ \n#$ -o "+tmpfoldname+"/outputs/\n#$ -q normal.q\n#$ -t 1-"+str(i-1)+"\n#$ -tc 4\n#$ -b y \n#$ -V \n#$ -S /bin/bash \n/home/galaxydev/galaxy/tools/SouthGreen/HMMer/hmmsearch "+profil_path+" ./base${SGE_TASK_ID} > ./outputs/output${SGE_TASK_ID}"
array_file=open(tmpfoldname+'/jobs.sge','w')
array_file.write(jobarray)
array_file.close()
log.close()

os.system("chmod +x "+tmpfoldname+"/jobs.sge")
os.system("qsub -sync y "+tmpfoldname+"/jobs.sge")
os.chdir(tmpfoldname)
os.system("tar -cf "+output_path+" outputs/output* 2> "+tmpfoldname+"/log.txt")   
os.system("tar -rf "+output_path+" bases.txt 2> "+tmpfoldname+"/log.txt")

os.chdir("../")
os.system("rm -rf "+tmpfoldname)
