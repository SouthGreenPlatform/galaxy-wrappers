#!/usr/bin/env python

"""
Input: fasta (input file), tabular (output file), int (truncation of id), int (columns from description)
Output: tabular
format convert: fasta to tabular
"""

import csv, re, sys, os, json
def stop_err( msg ):
    sys.stderr.write( msg )
    sys.exit()

def __main__():
    infile = sys.argv[1]
    outfile = sys.argv[2]
    home= os.path.expanduser('~')
    data_exp = home+"/galaxy/tool-data/SouthGreen/Coffee_genes_expression_all_organs_rpkm_hub.txt"
    fichier=open(infile,"r")
    destination = open(outfile,"w")
    data=open(data_exp,"r")
    COFCA_file=open("/bank/genfam/genome_data/COFCA/COFCA-GENOSCOPE1-sequence_feature-locus_tag-genfam.json","r")
    COFCA_dic=json.loads(COFCA_file.read())
    COFCA_file.close()
    
    data_dict=csv.DictReader(data, delimiter='\t')    
    gen_dict={}
    
    for line in data_dict:
        gen_dict[line["Gene"].upper()] = line
        
    
    for k in gen_dict:
        del gen_dict[k]["Gene"]
    
    list_name=[]
    for line in fichier :
        res = re.search(r"^>", line)
        if(res!=None):
            line.strip()
            if("COFCA" in line):
                name=re.sub(r"^>","",line)
                name=re.sub("\n","",name)
                origname=name
                name=re.sub(r"_[A-Z]{5}","",name)
                nametr=re.sub(r"C[PG]","CT",name)
                print name
                transc_name=COFCA_dic["COFCA"][nametr]["mrna_name"]
                transc_name=re.sub("_","",transc_name)
                list_name.append(name+";"+transc_name+";"+origname)
            
    organs=data_dict.fieldnames
    print(list_name[0:10])
    destination.write("\t".join(organs)+"\n")
    
    for gene in list_name :
        print gene
        name= gene.split(";")[0]
        mrna= gene.split(";")[1]
        origname= gene.split(";")[2]
        #name=re.sub(r"[pg]","t",gene)
        #name=re.sub(r"_","",name)
        #name=name+".1"
        try :
            ligne=[origname]
            for organe in organs[1:]:
                ligne.append(gen_dict[mrna][organe])
            destination.write("\t".join(ligne)+"\n")
        except:
            ligne=[origname]
            for y in range(1,len(organs)):
                ligne.append("UN")
            destination.write("\t".join(ligne)+"\n")
    destination.close()
    fichier.close()
    data.close()

if __name__ == "__main__" : __main__()
