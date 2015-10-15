#!/usr/bin/env python

"""
Input: fasta (input file), tabular (output file), int (truncation of id), int (columns from description)
Output: tabular
format convert: fasta to tabular
"""

import sys, os, re, subprocess
def stop_err( msg ):
    sys.stderr.write( msg )
    sys.exit()

#def __main__():
specie=sys.argv[1]
if(specie!="NONE"):
    test=os.popen("cat /SATA/galaxy/database/files/000/101/dataset_101131.dat | grep -w "+specie).read()
    if(test!=""):
        infile = sys.argv[2]
        outfile = sys.argv[3]
        expression=r"^>[^|\s]+[\s|]"
        fichier=open(infile,"r")
        destination = open(outfile,"w")
        for line in fichier :
            res = re.search(expression, line)
            if(res!=None):
                res2=res.group(0)[0:-1]
                formatted=re.search(r"_[A-Z]{5}$",res2)
                if(formatted!=None):
                    #new_line=re.sub(res2+" ",res2+"|",line)
                    new_line=res2+"\r\n"
                    new_line=re.sub(r"[ ,;]","_",new_line)
                    new_line=re.sub("__","_",new_line)
                    new_line=re.sub(r"[ _;| ,]*\r\n","\r\n",new_line)
                    destination.write(new_line)
                else:
                    #new_line=re.sub(res2+" ",res2+"_"+specie+"|",line)
                    new_line = res2+"_"+specie+"\r\n"
                    if(new_line==line):
                        new_line=re.sub(res2,res2+"_"+specie,line)
                    new_line=re.sub(r"[ ,;]","_",new_line)
                    new_line=re.sub("__","_",new_line)
                    new_line=re.sub(r"[ _;| ,]*\r\n","\r\n",new_line)
                    destination.write(new_line)
            else:
                destination.write(line)
        destination.close()
        fichier.close()
    else:
        stop_err("Incorrect Specie Code \n")
        
else:
    infile = sys.argv[2]
    outfile = sys.argv[3]
    # expression=r"^>[^|\s]+[\s|]"
    expression=r"\A>[0-z]*[^|\s\/]"
    fichier=open(infile,"r")
    destination = open(outfile,"w")
    
    for line in fichier :
		res = re.search(expression, line)
		if(res!=None):
			# name=res.group(0)[0:-1]
			name=res.group(0)
			if re.search(r">[A-Z]{5}_", line, flags=0) and re.search(r">GSMUA_", line, flags=0)==None and re.search(r">GRMZM_", line, flags=0)==None:
				spe_code=name[1:6]
				name=name[0]+name[7:]
				name=name+"_"+spe_code
			new_line=name+"\n"
			destination.write(new_line)
		else:
			destination.write(line)
    destination.close()
    fichier.close()

#if __name__ == "__main__" : __main__()
