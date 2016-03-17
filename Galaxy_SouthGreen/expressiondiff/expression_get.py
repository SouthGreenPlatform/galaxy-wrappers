#!/usr/bin/env python
# -*- coding: utf8 -*-

import sys, csv, re
import json

def import_exp_file(filename):
    file= open(filename, 'r')
    dic = csv.DictReader(file, delimiter='\t')
    dic_gene={}
    feature_name=filename[:-4]
    for row in dic:
        dic_gene[row["Gene.symbol"].upper()] = row
    return [feature_name, dic_gene]

input=sys.argv[1]
output= sys.argv[2]
fasta = open(input,'r')
res_file = open(output,'w')
rap_msu_file=open("/bank/genfam/RAP-MSU.txt",'r')
dic_rap_to_msu={}
dic_msu_to_rap={}

for line in rap_msu_file :
    rap_msu = line.split('\t')
    rap=rap_msu[0].strip().upper()
    msu=rap_msu[1].split(',')
    dic_rap_to_msu[rap]=msu
    for i in msu :
        dic_msu_to_rap[i.strip().upper()]=rap


rap_msu_file.close()


ARATH_file=open("/bank/genfam/genome_data/ARATH/ARATH-TAIR10-sequence_feature-locus_tag-genfam.json",'r')
ARATH_dic=json.loads(ARATH_file.read())
ARATH_file.close()

ORYSJ_file=open("/bank/genfam/genome_data/ORYSJ/ORYSJ-MSU7-sequence_feature-locus_tag-genfam.json",'r')
ORYSJ_dic=json.loads(ORYSJ_file.read())
ORYSJ_file.close()

[feat_ARATH_drough , ARATH_drough_genes] = import_exp_file("/bank/genfam/expression/ARATH/GEO_GSE11538_rosette_drough.txt")
[feat_ARATH_salt_WS , ARATH_salt_WS_genes] = import_exp_file("/bank/genfam/expression/ARATH/GSE16765 WS salt stress.txt")
[feat_ARATH_salt_col , ARATH_salt_col_genes] = import_exp_file("/bank/genfam/expression/ARATH/GSE16765 col salt stress.txt")
[feat_ORYSJ_salt_IR29 , ORYSJ_salt_IR29_genes] = import_exp_file("/bank/genfam/expression/ORYSJ/GEO2salinity stress response IR29.txt")
[feat_ORYSJ_salt_FL478 , ORYSJ_salt_FL478_genes] = import_exp_file("/bank/genfam/expression/ORYSJ/GEO2salinity stress response FL478.txt")

[feat_Drough_responsiveness_rice_all , Drough_responsiveness_rice_all_genes] = import_exp_file("/bank/genfam/expression/ORYSJ/GSE26280_Drough_responsiveness_rice_all_DK151.txt")
[feat_Drough_responsiveness_rice_leaves , Drough_responsiveness_rice_leaves_genes] = import_exp_file("/bank/genfam/expression/ORYSJ/GSE26280_Drough_responsiveness_rice_leaves_DK151.txt")
[feat_Drough_responsiveness_rice_roots , Drough_responsiveness_rice_roots_genes] = import_exp_file("/bank/genfam/expression/ORYSJ/GSE26280_Drough_responsiveness_rice_roots_DK151.txt")
[feat_Drough_responsiveness_rice_panicles, Drough_responsiveness_rice_panicles_genes] = import_exp_file("/bank/genfam/expression/ORYSJ/GSE26280_Drough_responsiveness_rice_panicles_DK151.txt")

res_file.write("Name\t annotation\t GSE11538_ARATH_drough\t GSE16765_ARATH_salt_col \t GSE16765_ARATH_salt_WS+ \t ORYSJ_salt_IR29 \t ORYSJ_salt_FL478_val \t GSE26280_Drough_responsiveness_rice_all \t GSE26280_Drough_responsiveness_rice_roots \t GSE26280_Drough_responsiveness_rice_leaves \t GSE26280_Drough_responsiveness_rice_panicles \n")

for line in fasta:
    if line[0]==">":
        prot_name=line.strip()[1:]
        specie_code=prot_name[-5:]
        locus_tag= prot_name[:-6]
        annot=[]
        if specie_code=="ARATH":
            gene_name=ARATH_dic["ARATH"][locus_tag]["gene_name"]
            gene_symbol=ARATH_dic["ARATH"][gene_name]["gene_symbol"]
            original_name=prot_name
            annot=[gene_symbol]
            #print gene_name
            if gene_name in ARATH_drough_genes.keys():
                ARATH_drough_val=ARATH_drough_genes[gene_name]["logFC"]
                if ARATH_drough_genes[gene_name]["Gene.title"] not in annot:
                    annot.append(ARATH_drough_genes[gene_name]["Gene.title"])
            elif  gene_symbol in ARATH_drough_genes.keys(): 
                ARATH_drough_val=ARATH_drough_genes[gene_symbol]["logFC"]
                if ARATH_drough_genes[gene_symbol]["Gene.title"] not in annot:
                    annot.append(ARATH_drough_genes[gene_symbol]["Gene.title"])
            else :
                ARATH_drough_val="NA"
            if gene_name in ARATH_salt_col_genes.keys() :
                ARATH_salt_col_val=ARATH_salt_col_genes[gene_name]["logFC"]
                if ARATH_salt_col_genes[gene_name]["Gene.title"] not in annot:
                    annot.append(ARATH_salt_col_genes[gene_name]["Gene.title"])
            elif  gene_symbol in ARATH_salt_col_genes.keys() : 
                ARATH_salt_col_val=ARATH_salt_col_genes[gene_symbol]["logFC"]
                if ARATH_salt_col_genes[gene_symbol]["Gene.title"] not in annot:
                    annot.append(ARATH_salt_col_genes[gene_symbol]["Gene.title"])
            else :
                ARATH_salt_col_val="NA"
            if gene_name in ARATH_salt_WS_genes.keys() :
                ARATH_salt_WS_val=ARATH_salt_WS_genes[gene_name]["logFC"]
                if ARATH_salt_WS_genes[gene_name]["Gene.title"] not in annot:
                    annot.append(ARATH_salt_WS_genes[gene_name]["Gene.title"])
            elif  gene_symbol in ARATH_salt_WS_genes.keys() : 
                ARATH_salt_WS_val=ARATH_salt_WS_genes[gene_symbol]["logFC"]
                if ARATH_salt_WS_genes[gene_symbol]["Gene.title"] not in annot:
                    annot.append(ARATH_salt_WS_genes[gene_symbol]["Gene.title"])
            else :
                ARATH_salt_WS_val="NA"
            ORYSJ_salt_IR29_val="NA"
            ORYSJ_salt_FL478_val="NA"
            Drough_responsiveness_rice_all_val="NA"
            Drough_responsiveness_rice_roots_val="NA"
            Drough_responsiveness_rice_leaves_val="NA"
            Drough_responsiveness_rice_panicles_val="NA"
            res_file.write(original_name+"\t"+'|'.join(annot)+"\t"+ ARATH_drough_val + "\t"+ARATH_salt_col_val+ "\t"+ARATH_salt_WS_val+ "\t"+ORYSJ_salt_IR29_val+"\t"+ORYSJ_salt_FL478_val+"\t"+ Drough_responsiveness_rice_all_val +"\t"+ Drough_responsiveness_rice_roots_val +"\t"+ Drough_responsiveness_rice_leaves_val +"\t"+ Drough_responsiveness_rice_panicles_val +"\n")
        elif specie_code=="ORYSJ" or specie_code=="ORYSA":
            print locus_tag
            gene_name=ORYSJ_dic["ORYSJ"][locus_tag]["gene_name"]
            print gene_name
            try:
                gene_symbol=ORYSJ_dic["ORYSJ"][gene_name]["gene_symbol"]
            except :
                gene_symbol=""
            original_name=prot_name
            msu_name=gene_name+".1"
            annot.append(locus_tag)
            annot.append(gene_name)
            #if ():      
            try:
                gene_name= dic_msu_to_rap[msu_name]
            except:
                gene_name=msu_name
            annot=[gene_symbol]
            #print gene_name
            ARATH_drough_val="NA"
            ARATH_salt_col_val="NA"
            ARATH_salt_WS_val="NA"
            if gene_name in ORYSJ_salt_IR29_genes.keys():
                ORYSJ_salt_IR29_val=ORYSJ_salt_IR29_genes[gene_name]["logFC"]
                if ORYSJ_salt_IR29_genes[gene_name]["Gene.title"] not in annot:
                    annot.append(ORYSJ_salt_IR29_genes[gene_name]["Gene.title"])
            elif gene_symbol in ORYSJ_salt_IR29_genes.keys() : 
                ORYSJ_salt_IR29_val=ORYSJ_salt_IR29_genes[gene_symbol]["logFC"]
                if ORYSJ_salt_IR29_genes[gene_symbol]["Gene.title"] not in annot:
                    annot.append(ORYSJ_salt_IR29_genes[gene_symbol]["Gene.title"])
            else :
                ORYSJ_salt_IR29_val="NA"
            if gene_name in ORYSJ_salt_FL478_genes.keys() :
                ORYSJ_salt_FL478_val=ORYSJ_salt_FL478_genes[gene_name]["logFC"]
                if ORYSJ_salt_FL478_genes[gene_name]["Gene.title"] not in annot:
                    annot.append(ORYSJ_salt_FL478_genes[gene_name]["Gene.title"])
            elif gene_symbol in ORYSJ_salt_FL478_genes.keys() : 
                ORYSJ_salt_FL478_val=ORYSJ_salt_FL478_genes[gene_symbol]["logFC"]
                if ORYSJ_salt_FL478_genes[gene_symbol]["Gene.title"] not in annot:
                    annot.append(ORYSJ_salt_FL478_genes[gene_symbol]["Gene.title"])
            else:
                ORYSJ_salt_FL478_val="NA"
            if gene_name in Drough_responsiveness_rice_all_genes.keys() :
                Drough_responsiveness_rice_all_val=Drough_responsiveness_rice_all_genes[gene_name]["logFC"]
                if Drough_responsiveness_rice_all_genes[gene_name]["Gene.title"] not in annot:
                    annot.append(Drough_responsiveness_rice_all_genes[gene_name]["Gene.title"])
            elif gene_symbol in Drough_responsiveness_rice_all_genes.keys() :
                Drough_responsiveness_rice_all_val=Drough_responsiveness_rice_all_genes[gene_symbol]["logFC"]
                if Drough_responsiveness_rice_all_genes[gene_symbol]["Gene.title"] not in annot:
                    annot.append(Drough_responsiveness_rice_all_genes[gene_symbol]["Gene.title"])
            else :
                Drough_responsiveness_rice_all_val="NA"
            if gene_name in Drough_responsiveness_rice_roots_genes.keys()  :
                Drough_responsiveness_rice_roots_val=Drough_responsiveness_rice_roots_genes[gene_name]["logFC"]
                if Drough_responsiveness_rice_roots_genes[gene_name]["Gene.title"] not in annot:
                    annot.append(Drough_responsiveness_rice_roots_genes[gene_name]["Gene.title"])
            elif gene_symbol in Drough_responsiveness_rice_roots_genes.keys() : 
                Drough_responsiveness_rice_roots_val=Drough_responsiveness_rice_roots_genes[gene_symbol]["logFC"]
                if Drough_responsiveness_rice_roots_genes[gene_symbol]["Gene.title"] not in annot:
                    annot.append(Drough_responsiveness_rice_roots_genes[gene_symbol]["Gene.title"])
            else : 
                Drough_responsiveness_rice_roots_val="NA"
            if gene_name in Drough_responsiveness_rice_leaves_genes.keys() :
                Drough_responsiveness_rice_leaves_val=Drough_responsiveness_rice_leaves_genes[gene_name]["logFC"]
                if Drough_responsiveness_rice_leaves_genes[gene_name]["Gene.title"] not in annot:
                    annot.append(Drough_responsiveness_rice_leaves_genes[gene_name]["Gene.title"])
            elif gene_symbol in Drough_responsiveness_rice_leaves_genes.keys() : 
                Drough_responsiveness_rice_leaves_val=Drough_responsiveness_rice_leaves_genes[gene_symbol]["logFC"]
                if Drough_responsiveness_rice_leaves_genes[gene_symbol]["Gene.title"] not in annot:
                    annot.append(Drough_responsiveness_rice_leaves_genes[gene_symbol]["Gene.title"])
            else : 
                Drough_responsiveness_rice_leaves_val="NA"
            if gene_name in Drough_responsiveness_rice_panicles_genes.keys():
                Drough_responsiveness_rice_panicles_val=Drough_responsiveness_rice_panicles_genes[gene_name]["logFC"]
                if Drough_responsiveness_rice_panicles_genes[gene_name]["Gene.title"] not in annot:
                    annot.append(Drough_responsiveness_rice_panicles_genes[gene_name]["Gene.title"])
            elif gene_symbol in Drough_responsiveness_rice_panicles_genes.keys() : 
                Drough_responsiveness_rice_panicles_val=Drough_responsiveness_rice_panicles_genes[gene_symbol]["logFC"]
                if Drough_responsiveness_rice_panicles_genes[gene_symbol]["Gene.title"] not in annot:
                    annot.append(Drough_responsiveness_rice_panicles_genes[gene_symbol]["Gene.title"])
            else : 
                Drough_responsiveness_rice_panicles_val="NA"
            res_file.write(original_name+"\t"+'|'.join(annot)+"\t"+ ARATH_drough_val + "\t"+ARATH_salt_col_val+ "\t"+ARATH_salt_WS_val+ "\t"+ORYSJ_salt_IR29_val+"\t"+ORYSJ_salt_FL478_val+"\t"+ Drough_responsiveness_rice_all_val +"\t"+ Drough_responsiveness_rice_roots_val +"\t"+ Drough_responsiveness_rice_leaves_val +"\t"+ Drough_responsiveness_rice_panicles_val +"\n")
        #print specie_code
        #if specie_code=="ARATH":
            

fasta.close()
res_file.close()
