#!/usr/bin/env python

"""
Input: fasta (input file), tabular (output file), int (truncation of id), int (columns from description)
Output: tabular
format convert: fasta to tabular
"""

import csv, re, sys
def stop_err( msg ):
    sys.stderr.write( msg )
    sys.exit()

def __main__():
    infile = sys.argv[1]
    outfile = sys.argv[2]
    data_exp = " ~/galaxy/tool-data/SouthGreen/Coffee_genes_expression_all_organs_rpkm_hub.txt"
    fichier=open(infile,"r")
    destination = open(outfile,"w")
    data=open(data_exp,"r")
    
    data_dict=csv.DictReader(data, delimiter='\t')    
    gen_dict={}
    
    for line in data_dict:
        gen_dict[line["Gene"]] = line
        
    
    for k in gen_dict:
        del gen_dict[k]["Gene"]
    
    list_name=[]
    for line in fichier :
        res = re.search(r"^>", line)
        if(res!=None):
            line.strip()
            if("COFCA" in line):
                name=re.sub(r"_[A-Z]{5}","",line)
                name=re.sub(r"^>","",name)
                name=re.sub("\n","",name)
                print name
                list_name.append(name)
            
    organs=data_dict.fieldnames
    print(list_name[0:10])
    destination.write(",".join(organs)+"\n")
    
    for gene in list_name :
        print gene
        name=re.sub(r"[pg]","t",gene)
        #name=re.sub(r"_","",name)
        #name=name+".1"
        try :
            ligne=[gene]
            for organe in organs[1:]:
                ligne.append(gen_dict[name][organe])
            destination.write(",".join(ligne)+"\n")
        except:
            ligne=[gene]
            for y in range(1,len(organs)):
                ligne.append("UN")
            destination.write(",".join(ligne)+"\n")
    destination.close()
    fichier.close()
    data.close()

if __name__ == "__main__" : __main__()
