#!/usr/bin/env python
# -*- coding: utf8 -*-


from multiprocessing import Process, Manager
import os, sys, random, datetime,re
import json

fasta_path = sys.argv[1]
output_path=sys.argv[2]


fasta_file=open(fasta_path,'r')
output_file=open(output_path,'w')
# rap_msu_file=open("/bank/genfam/RAP-MSU.txt",'r')
# dic_rap_to_msu={}
# dic_msu_to_rap={}

output_file.write("Protein_id\tProtein_name\tGene_Name\tGene_symbol\tUniprot_entry\tUniprot_product\n")
dic_of_dic={}

ARATH_file=open("/bank/genfam/genome_data/ARATH/ARATH-TAIR10-sequence_feature-locus_tag-genfam.json",'r')
ARATH_dic=json.loads(ARATH_file.read())
dic_of_dic["ARATH"]=ARATH_dic
ARATH_file.close()

BRADI_file=open("/bank/genfam/genome_data/BRADI/BRADI-JGI1-sequence_feature-locus_tag-genfam.json",'r')
BRADI_dic=json.loads(BRADI_file.read())
dic_of_dic["BRADI"]=BRADI_dic
BRADI_file.close()

COFCA_file=open("/bank/genfam/genome_data/COFCA/COFCA-GENOSCOPE1-sequence_feature-locus_tag-genfam.json",'r')
COFCA_dic=json.loads(COFCA_file.read())
dic_of_dic["COFCA"]=COFCA_dic
COFCA_file.close()

ELAGV_file=open("/bank/genfam/genome_data/ELAGV/ELAGV-MPOB5-sequence_feature-locus_tag-genfam.json",'r')
ELAGV_dic=json.loads(ELAGV_file.read())
dic_of_dic["ELAGV"]=ELAGV_dic
ELAGV_file.close()

ELAGV_REF_file=open("/bank/genfam/genome_data/REF_ELAGV/ELAGV-NCBI100-sequence_feature-locus_tag-genfam.json",'r')
ELAGV_REF_dic=json.loads(ELAGV_REF_file.read())
dic_of_dic["REF_ELAGV"]=ELAGV_REF_dic
ELAGV_REF_file.close()

EUCGR_REF_file=open("/bank/genfam/genome_data/REF_EUCGR/EUCGR-NCBI100-sequence_feature-locus_tag-genfam.json",'r')
EUCGR_REF_dic=json.loads(EUCGR_REF_file.read())
dic_of_dic["REF_EUCGR"]=EUCGR_REF_dic
EUCGR_REF_file.close()

GLYMA_file=open("/bank/genfam/genome_data/GLYMA/GLYMA-JGI1-sequence_feature-locus_tag-genfam.json",'r')
GLYMA_dic=json.loads(GLYMA_file.read())
dic_of_dic["GLYMA"]=GLYMA_dic
GLYMA_file.close()

GOSRA_file=open("/bank/genfam/genome_data/GOSRA/GOSRA-JGI1-sequence_feature-locus_tag-genfam.json",'r')
GOSRA_dic=json.loads(GOSRA_file.read())
dic_of_dic["GOSRA"]=GOSRA_dic
GOSRA_file.close()

LOTJA_file=open("/bank/genfam/genome_data/LOTJA/LOTJA-Kazusa2.5-sequence_feature-locus_tag-genfam.json",'r')
LOTJA_dic=json.loads(LOTJA_file.read())
dic_of_dic["LOTJA"]=LOTJA_dic
LOTJA_file.close()

MAIZE_file=open("/bank/genfam/genome_data/MAIZE/MAIZE-MGDB5b60-sequence_feature-locus_tag-genfam.json",'r')
MAIZE_dic=json.loads(MAIZE_file.read())
dic_of_dic["MAIZE"]=MAIZE_dic
MAIZE_file.close()

MALDO_file=open("/bank/genfam/genome_data/MALDO/MALDO-JGI1-sequence_feature-locus_tag-genfam.json",'r')
MALDO_dic=json.loads(MALDO_file.read())
dic_of_dic["MALDO"]=MALDO_dic
MALDO_file.close()

MANES_file=open("/bank/genfam/genome_data/MANES/MANES-JGI4.1-sequence_feature-locus_tag-genfam.json",'r')
MANES_dic=json.loads(MANES_file.read())
dic_of_dic["MANES"]=MANES_dic
MANES_file.close()

MEDTR_file=open("/bank/genfam/genome_data/MEDTR/MEDTR-Mt3.5v5-sequence_feature-locus_tag-genfam.json",'r')
MEDTR_dic=json.loads(MEDTR_file.read())
dic_of_dic["MEDTR"]=MEDTR_dic
MEDTR_file.close()

MUSAC_file=open("/bank/genfam/genome_data/MUSAC/MUSAC-Musa1-sequence_feature-locus_tag-genfam.json",'r')
MUSAC_dic=json.loads(MUSAC_file.read())
dic_of_dic["MUSAC"]=MUSAC_dic
MUSAC_file.close()

ORYSI_file=open("/bank/genfam/genome_data/ORYSI/ORYSI-BGI2-sequence_feature-locus_tag-genfam.json",'r')
ORYSI_dic=json.loads(ORYSI_file.read())
dic_of_dic["ORYSI"]=ORYSI_dic
ORYSI_file.close()

ORYSJ_file=open("/bank/genfam/genome_data/ORYSJ/ORYSJ-MSU7-sequence_feature-locus_tag-genfam.json",'r')
ORYSJ_dic=json.loads(ORYSJ_file.read())
dic_of_dic["ORYSJ"]=ORYSJ_dic
ORYSJ_file.close()

PHODC_file=open("/bank/genfam/genome_data/PHODC/PHODC-WCMCQ3-sequence_feature-locus_tag-genfam.json",'r')
PHODC_dic=json.loads(PHODC_file.read())
dic_of_dic["PHODC"]=PHODC_dic
PHODC_file.close()

POPTR_file=open("/bank/genfam/genome_data/POPTR/POPTR-JGI2-sequence_feature-locus_tag-genfam.json",'r')
POPTR_dic=json.loads(POPTR_file.read())
dic_of_dic["POPTR"]=POPTR_dic
POPTR_file.close()

RICCO_file=open("/bank/genfam/genome_data/RICCO/RICCO-JGI0.1-sequence_feature-locus_tag-genfam.json",'r')
RICCO_dic=json.loads(RICCO_file.read())
dic_of_dic["RICCO"]=RICCO_dic
RICCO_file.close()

SETIT_file=open("/bank/genfam/genome_data/SETIT/SETIT-JGI1-sequence_feature-locus_tag-genfam.json",'r')
SETIT_dic=json.loads(SETIT_file.read())
dic_of_dic["SETIT"]=SETIT_dic
SETIT_file.close()

SOLLC_file=open("/bank/genfam/genome_data/SOLLC/SOLLC-ITAG2.40-sequence_feature-locus_tag-genfam.json",'r')
SOLLC_dic=json.loads(SOLLC_file.read())
dic_of_dic["SOLLC"]=SOLLC_dic
SOLLC_file.close()

SOLTU_file=open("/bank/genfam/genome_data/SOLTU/SOLTU-JGI3-sequence_feature-locus_tag-genfam.json",'r')
SOLTU_dic=json.loads(SOLTU_file.read())
dic_of_dic["SOLTU"]=SOLTU_dic
SOLTU_file.close()

SORBI_file=open("/bank/genfam/genome_data/SORBI/SORBI-JGI1.4-sequence_feature-locus_tag-genfam.json",'r')
SORBI_dic=json.loads(SORBI_file.read())
dic_of_dic["SORBI"]=SORBI_dic
SORBI_file.close()

THECC_file=open("/bank/genfam/genome_data/THECC/THECC-COCOA1-sequence_feature-locus_tag-genfam.json",'r')
THECC_dic=json.loads(THECC_file.read())
dic_of_dic["THECC"]=THECC_dic
THECC_file.close()


VITVI_file=open("/bank/genfam/genome_data/VITVI/VITVI-GENOSCOPE1-sequence_feature-locus_tag-genfam.json",'r')
VITVI_dic=json.loads(VITVI_file.read())
dic_of_dic["VITVI"]=VITVI_dic
VITVI_file.close()


fasta=fasta_file.readlines()
for line in fasta:
	if line[0]==">":
	    prot_id=line[1:].strip()
	    print prot_id
	    specode=prot_id.split("_")[-1]
	    if re.search("REF_",prot_id):
	        specode="REF_"+specode
	    print specode
	    locus_tag="_".join(prot_id.split("_")[:-1])
	    print locus_tag
	    if "gene_name" in dic_of_dic[specode][specode][locus_tag].keys():
	        gene_name=dic_of_dic[specode][specode][locus_tag]["gene_name"]
	    else : 
	        gene_name=""
	    if "polypeptide_name" in dic_of_dic[specode][specode][locus_tag].keys():
	        prot_name=dic_of_dic[specode][specode][locus_tag]["polypeptide_name"]
	    else : 
	        prot_name=""
	    if "gene_symbol" in dic_of_dic[specode][specode][gene_name].keys():
	        gene_symbol=dic_of_dic[specode][specode][gene_name]["gene_symbol"]
	    else : 
	        gene_symbol=""
	    if "uniprot_entry" in dic_of_dic[specode][specode][gene_name].keys():
	        uniprot_entry=dic_of_dic[specode][specode][gene_name]["uniprot_entry"]
	    else : 
	        uniprot_entry=""
	    if "uniprot_product" in dic_of_dic[specode][specode][gene_name].keys():
	        uniprot_product=dic_of_dic[specode][specode][gene_name]["uniprot_product"]
	    else : 
	        uniprot_product=""
	    output_file.write(prot_id+"\t"+prot_name+"\t"+gene_name+"\t"+gene_symbol+"\t"+uniprot_entry+"\t"+uniprot_product+"\n")
	    
fasta_file.close()
output_file.close()