#!/usr/bin/env python
# -*- coding: utf8 -*-

import sys, csv, re


def import_exp_file(filename):
    file= open(filename, 'r')
    dic = csv.DictReader(file, delimiter='\t')
    dic_gene={}
    feature_name=filename[:-4]
    for row in dic:
        dic_gene[row["Gene.symbol"]] = row
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
    rap=rap_msu[0].strip()
    msu=rap_msu[1].split(',')
    dic_rap_to_msu[rap]=msu
    for i in msu :
        dic_msu_to_rap[i.strip()]=rap

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
        gene_name=line.strip()[1:]
        specie_code=gene_name[-5:]
        locus_tag= gene_name[:-6]
        annot=[]
        if specie_code=="ARATH":
            original_name=gene_name
            gene_name=gene_name.split('.')[0]
            #print gene_name
            try :
                ARATH_drough_val=ARATH_drough_genes[gene_name]["logFC"]
                if ARATH_drough_genes[gene_name]["Gene.title"] not in annot:
                    annot.append(ARATH_drough_genes[gene_name]["Gene.title"])
            except : 
                ARATH_drough_val="NA"
            try :
                ARATH_salt_col_val=ARATH_salt_col_genes[gene_name]["logFC"]
                if ARATH_salt_col_genes[gene_name]["Gene.title"] not in annot:
                    annot.append(ARATH_salt_col_genes[gene_name]["Gene.title"])
            except : 
                ARATH_salt_col_val="NA"
            try :
                ARATH_salt_WS_val=ARATH_salt_WS_genes[gene_name]["logFC"]
                if ARATH_salt_WS_genes[gene_name]["Gene.title"] not in annot:
                    annot.append(ARATH_salt_WS_genes[gene_name]["Gene.title"])
            except : 
                ARATH_salt_WS_val="NA"
            ORYSJ_salt_IR29_val="NA"
            ORYSJ_salt_FL478_val="NA"
            Drough_responsiveness_rice_all_val="NA"
            Drough_responsiveness_rice_roots_val="NA"
            Drough_responsiveness_rice_leaves_val="NA"
            Drough_responsiveness_rice_panicles_val="NA"
            res_file.write(original_name+"\t"+'|'.join(annot)+"\t"+ ARATH_drough_val + "\t"+ARATH_salt_col_val+ "\t"+ARATH_salt_WS_val+ "\t"+ORYSJ_salt_IR29_val+"\t"+ORYSJ_salt_FL478_val+"\t"+ Drough_responsiveness_rice_all_val +"\t"+ Drough_responsiveness_rice_roots_val +"\t"+ Drough_responsiveness_rice_leaves_val +"\t"+ Drough_responsiveness_rice_panicles_val +"\n")
        elif specie_code=="ORYSJ" or specie_code=="ORYSA":
            msu_name=locus_tag
            original_name=gene_name
            #if ():
                
            try:
                gene_name= dic_msu_to_rap[msu_name]
            except:
                gene_name=msu_name
            annot=[gene_name]
            #print gene_name
            ARATH_drough_val="NA"
            ARATH_salt_col_val="NA"
            ARATH_salt_WS_val="NA"
            try :
                ORYSJ_salt_IR29_val=ORYSJ_salt_IR29_genes[gene_name]["logFC"]
                if ORYSJ_salt_IR29_genes[gene_name]["Gene.title"] not in annot:
                    annot.append(ORYSJ_salt_IR29_genes[gene_name]["Gene.title"])
            except : 
                ORYSJ_salt_IR29_val="NA"
            try :
                ORYSJ_salt_FL478_val=ORYSJ_salt_FL478_genes[gene_name]["logFC"]
                if ORYSJ_salt_FL478_genes[gene_name]["Gene.title"] not in annot:
                    annot.append(ORYSJ_salt_FL478_genes[gene_name]["Gene.title"])
            except : 
                ORYSJ_salt_FL478_val="NA"
            try :
                Drough_responsiveness_rice_all_val=Drough_responsiveness_rice_all_genes[gene_name]["logFC"]
                if Drough_responsiveness_rice_all_genes[gene_name]["Gene.title"] not in annot:
                    annot.append(Drough_responsiveness_rice_all_genes[gene_name]["Gene.title"])
            except : 
                Drough_responsiveness_rice_all_val="NA"
            try :
                Drough_responsiveness_rice_roots_val=Drough_responsiveness_rice_roots_genes[gene_name]["logFC"]
                if Drough_responsiveness_rice_roots_genes[gene_name]["Gene.title"] not in annot:
                    annot.append(Drough_responsiveness_rice_roots_genes[gene_name]["Gene.title"])
            except : 
                Drough_responsiveness_rice_roots_val="NA"
            try :
                Drough_responsiveness_rice_leaves_val=Drough_responsiveness_rice_leaves_genes[gene_name]["logFC"]
                if Drough_responsiveness_rice_leaves_genes[gene_name]["Gene.title"] not in annot:
                    annot.append(Drough_responsiveness_rice_leaves_genes[gene_name]["Gene.title"])
            except : 
                Drough_responsiveness_rice_leaves_val="NA"
            try :
                Drough_responsiveness_rice_panicles_val=Drough_responsiveness_rice_panicles_genes[gene_name]["logFC"]
                if Drough_responsiveness_rice_panicles_genes[gene_name]["Gene.title"] not in annot:
                    annot.append(Drough_responsiveness_rice_panicles_genes[gene_name]["Gene.title"])
            except : 
                Drough_responsiveness_rice_panicles_val="NA"
            res_file.write(original_name+"\t"+'|'.join(annot)+"\t"+ ARATH_drough_val + "\t"+ARATH_salt_col_val+ "\t"+ARATH_salt_WS_val+ "\t"+ORYSJ_salt_IR29_val+"\t"+ORYSJ_salt_FL478_val+"\t"+ Drough_responsiveness_rice_all_val +"\t"+ Drough_responsiveness_rice_roots_val +"\t"+ Drough_responsiveness_rice_leaves_val +"\t"+ Drough_responsiveness_rice_panicles_val +"\n")
        #print specie_code
        #if specie_code=="ARATH":
            

fasta.close()
res_file.close()
