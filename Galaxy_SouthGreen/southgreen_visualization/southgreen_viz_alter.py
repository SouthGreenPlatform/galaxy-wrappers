#!/usr/bin/env python
# -*- coding: utf8 -*-

from multiprocessing import Process, Manager
import os, sys, random, datetime


tree_path = sys.argv[1]
al_path=sys.argv[2]
ideven_path=sys.argv[3]
exp_path=sys.argv[4]
genfam_path=sys.argv[5]

curdir=os.getcwd()
random.seed()
numb=str(random.randint(1,1000000))
tmpfoldname=curdir+"/"+numb
os.system("mkdir "+tmpfoldname)
os.system("chmod 777 "+tmpfoldname)

os.system("tar -c "+genfam_path+ " 2> /dev/null")

if tree_path!="notree":
    os.system("cp "+tree_path+" "+tmpfoldname+"/tree_file.txt")
    os.system("tar -rvf "+genfam_path+" -C "+tmpfoldname+"/ tree_file.txt 2> /dev/null ")
if al_path!="noal":
    os.system("cp "+al_path+" "+tmpfoldname+"/al_file.txt")
    os.system("tar -rvf "+genfam_path+" -C "+tmpfoldname+"/ al_file.txt 2> /dev/null ")
if ideven_path!="noideven":
    os.system("cp "+ideven_path+" "+tmpfoldname+"/ideven_file.txt")
    os.system("tar -rvf "+genfam_path+" -C "+tmpfoldname+"/ ideven_file.txt 2> /dev/null ")
if exp_path!="noexp":
    os.system("cp "+exp_path+" "+tmpfoldname+"/exp_file.txt")
    os.system("tar -rvf "+genfam_path+" -C "+tmpfoldname+"/ exp_file.txt 2> /dev/null ")

# os.system("tar -cvf "+genfam_path+" "+tmpfoldname+"*")
os.system("rm -rf "+tmpfoldname)