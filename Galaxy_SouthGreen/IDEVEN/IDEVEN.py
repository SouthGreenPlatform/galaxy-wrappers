#!/usr/bin/env python
# -*- coding: utf8 -*-

from multiprocessing import Process, Manager
import os, sys, random, datetime
import json
import pickle

family_path = sys.argv[1]
output_path=sys.argv[2]

# with open('/home/galaxydev/SouthGreen_tools/Galaxy_SouthGreen/IDEVEN/dictionnaire_alias.pickle', 'rb') as fichier:
    # mon_depickler = pickle.Unpickler(fichier)
    # dico_recupere = mon_depickler.load()

fam_file=open(family_path,'r')
fam=fam_file.readlines()
fam_file_out=open(family_path,'r')
speList=[]
for line in fam:
    if line[0]==">":
        if "REF" in line:
            species="REF_"+line.strip()[-5:]
        else:
            species=line.strip()[-5:]
        oldname=line.strip()[:-6]
        try: 
            newname=dico_recupere[line.strip()[-5:]][oldname]["polypeptide_name"]
            line=line.replace(oldname,newname)
        except: 
            line=line
        if species  not in speList : 
            speList.append(species)
            print species 
		
# print speList
fam_file.close()
fam_file_out.close()
curdir=os.getcwd()
random.seed()
numb=str(random.randint(1,1000000))
tmpfoldname=curdir+"/"+numb
os.system("mkdir "+tmpfoldname)
os.system("chmod 777 "+tmpfoldname)
i=0

directory=os.getcwd()
print directory

while i<len(speList) :
    j=i
    while j<len(speList) : 
#         log_file=open("/home/galaxydev/galaxy/tools/SouthGreen/IDEVEN/"+speList[i]+"_"+speList[j]+"_error.log",'w')
#        print ("qsub -V -q web.q -o  "+tmpfoldname+"/ouput.log -sync y -N IDEVEN_"+speList[i]+"_"+speList[j]+" -b y '/home/galaxydev/galaxy/tools/SouthGreen/IDEVEN/IDEVEN.sh "+speList[i]+" "+speList[j]+" "+family_path+" "+tmpfoldname+"/IDEVEN_"+speList[i]+"_"+speList[j]+".txt'")
        # print "qsub -V -q web.q -l mem_free=12G -o "+tmpfoldname+"/ouput.log -e "+tmpfoldname+"/error.log -sync y -N IDEVEN_"+speList[i]+"_"+speList[j]+" -b y '/home/galaxydev/galaxy/tools/SouthGreen/IDEVEN/IDEVEN.sh "+speList[i]+" "+speList[j]+" "+family_path+" "+tmpfoldname+"/IDEVEN_"+speList[i]+"_"+speList[j]+".txt' 2> "+tmpfoldname+"/"+speList[i]+"_"+speList[j]+"_error.log "
        if speList[i]!= speList[j]:
			print("qsub -V -q web.q -l mem_free=12G -o "+tmpfoldname+"/"+speList[i]+"_"+speList[j]+"_ouput.log -e "+tmpfoldname+"/"+speList[i]+"_"+speList[j]+"_error.log -sync y -N IDEVEN_"+speList[i]+"_"+speList[j]+" -b y '/home/galaxydev/galaxy/tools/SouthGreen/IDEVEN/IDEVEN.sh "+speList[i]+" "+speList[j]+" "+family_path+" "+tmpfoldname+"/IDEVEN_"+speList[i]+"_"+speList[j]+".txt' 2> "+tmpfoldname+"/"+speList[i]+"_"+speList[j]+"_error.log ")
			os.system("qsub -V -q web.q -l mem_free=12G -o "+tmpfoldname+"/"+speList[i]+"_"+speList[j]+"_ouput.log -e "+tmpfoldname+"/"+speList[i]+"_"+speList[j]+"_error.log -sync y -N IDEVEN_"+speList[i]+"_"+speList[j]+" -b y '/home/galaxydev/galaxy/tools/SouthGreen/IDEVEN/IDEVEN.sh "+speList[i]+" "+speList[j]+" "+family_path+" "+tmpfoldname+"/IDEVEN_"+speList[i]+"_"+speList[j]+".txt' 2> "+tmpfoldname+"/"+speList[i]+"_"+speList[j]+"_error.log ")
        j=j+1
    i=i+1

out_file=open(output_path,'w')
out_file.write("#Nom gene1\tNom Gene2\tevent\tKs\tMean Ks\tsize bloc\n")
out_file.close()
os.system("tail -q -n +2 "+tmpfoldname+"/IDEVEN_*.txt >> "+output_path)
os.system("rm -rf "+tmpfoldname)

