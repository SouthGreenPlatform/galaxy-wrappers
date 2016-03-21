#!/usr/bin/env python
# -*- coding: utf8 -*-

from multiprocessing import Process, Manager
import os, sys, random, datetime,re
import json

# Pour chaque organisme : regex générale du type de nomenclature, puis à chaque ligne du fichier, test si on trouve  une des correspondances, et selon les cas, ajout du code espèce, remplacement P/T/G et/ou transcrit 1 avec ".1" ou juste sans extension

remove_tr= sys.argv[1] # yes or no
fasta_path = sys.argv[2] # fasta a trier
output_path=sys.argv[3] #  fasta resultat

fasta_file=open(fasta_path,'r')
output_file=open(output_path,'w')
rap_msu_file=open("/bank/genfam/RAP-MSU.txt",'r') #fichier d'équivalence RAP/MSU pour le riz
dic_rap_to_msu={}
dic_msu_to_rap={}
# création des dictionnaires de traduction RAP vers MSU et MSU vers RAP
for line in rap_msu_file :
    rap_msu = line.split('\t')
    rap=rap_msu[0].strip()
    msu=rap_msu[1].split(',')
    dic_rap_to_msu[rap]=msu[0]
    for i in msu :
        dic_msu_to_rap[i.strip()]=rap

		
fasta=fasta_file.readlines()
# pat_ARATH_transc=r'(AT[0-9]+G[0-9]{5}\.[0-9]+[^_])'
# pat_ARATH=r'(AT[0-9]+G[0-9]{5}[^_\.])'
# pat_bradi_transc=r'(Bradi[0-9]+g[0-9]{5}\.[0-9]+[^_])'
# pat_bradi=r'(Bradi[0-9]+g[0-9]{5}[^_\.])'
# pat_glyma=r'(Glyma[0-9]+g[0-9]{5}[^_\.])'
# pat_glyma_transc=r'(Glyma[0-9]+g[0-9]{5}\.[0-9]+[^_])'
# pat_gosra=r'(Gorai\.[0-9]+G[0-9]{6}\.[0-9]+[^_])'
# pat_lotja=r'(LjSGA_[0-9]{6}\.[0-9]+[^_])' #pb id avec coge
# pat_maize=r'(GRMZM[0-9]+[G][0-9]{6}\_[TGP][0-9]+[^_])'
# pat_MALDO=r'(MD[CP][0-9]{10}(\.)*[0-9]*[^_])'
# pat_manes=r'(cassava4\.1_[0-9]{6}m[^_])'
# pat_medtr=r'(Medtr[0-9]+g[0-9]{6}[^_\.])'
# pat_medtr_transcr=r'(Medtr[0-9]+g[0-9]{6}[^_][0-9]+)'
# pat_musa=r'(GSMUA_Achr([0-9]+|(Un_random))[TGP][0-9]{5}\_[0-9]{3}[^_])'
# pat_orysj_MSU=r'((LOC_)?Os[0-9]{2}g[0-9]{5}(\.[0-9]+)?[^0-9])'
# pat_orysj_RAP=r'(Os[0-9]{2}g[0-9]{7}[a|b]?)'
# pat_ricco=r'([0-9]{5}.m[0-9]{6}(\.[0-9]+)*[^_])'
# pat_sollc=r'(Solyc[0-9]{2}g[0-9]{6}\.[0-9]+\.[0-9]+[^_])'
# pat_soltu=r'(PGSC[0-9]{4}DM[GTP][0-9]{9}(\.[0-9]+)*[^_])'
# pat_sorbi=r'(Sb[0-9]{2}g[0-9]{6}\.[0-9]+[^_])'
# pat_cucsa=r'(Cucsa\.[0-9]{6}(\.[0-9]+)?)'
# pat_citsi=r'(orange[0-9]+\.[0-9]+g[0-9]{6}m[^_])'
# pat_horvu=r'((CDS:)*MLOC_[0-9]{5}\.[0-9]+[^_])'
# pat_CAJCA=r'(C.cajan_[0-9]{5}(\[.[0-9]+)[^_])'
# pat_MUSBA=r'(ITC[0-9]{4}_Bchr[0-9]+_[TPG][0-9]+[^_])'
# pat_POPTR=r'(Potri\.[0-9]{3}G[0-9]{6}\.[1-9]+[^_])' #pb id avec coge
# pat_PHODC=r'(PDK_[0-9]{2}s[0-9]+[Lg][0-9]{3}[^_])'
# pat_thecc=r'(Tc[0-9]{2}_[tg][0-9]{6}[^_])' #pb id avec coge
# pat_theccbis=r'(TCM\_[0-9]{6}[^_])' #pb id avec coge
# pat_vitvi=r'(GSVIVT[0-9]{11}[^_])'
# pat_vitvibis=r'(LOC[0-9]{9}[^_])'

# gfpat_ARATH= "AT{1-9}G04860.1-PROTEIN"
# gfpat_COFCA=
# gfpat_BRADI=
# gfpat_ELAGV=
# gfpat_ELAGV_REF=
# gfpat_EUCGR_REF=
# gfpat_GLYMA=
# gfpat_GOSRA=
# gfpat_LOTJA=
# gfpat_MAIZE=
# gfpat_MALDO=
# gfpat_MANES=
# gfpat_MEDTR=
# gfpat_MUSAC=
# gfpat_ORYSI=
# gfpat_ORYSJ=
# gfpat_RICCO=
# gfpat_PHODC=
# gfpat_POPTR=
# gfpat_SOLTU=
# gfpat_SETIT=
# gfpat_SOLLC=
# gfpat_VITVI=
# gfpat_SORBI=
# gfpat_THECC=
# 
# 
# # pour emélioration : mettre dico en mempoire vive, ou préprer des fichiers binaires des dico pour un acces plus rapide
# 
rap_msu_file=open("/bank/genfam/RAP-MSU.txt",'r') #fichier d'équivalence RAP/MSU pour le riz
dic_rap_to_msu={}
dic_msu_to_rap={}

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

fasta_file.close()


fasta_file=open(fasta_path,'r')
seqfasta=""
for line in fasta:
	if line[0]==">":
		name = line[1:].upper()
		# Si l'option de retrait des transcrits alternatifs est "yes" : on teste si la séquence fasta précédente est un transcrit alternatif, et dans ce cas on ne l'écrit pas dans le fichier de sortie
		#test
		if remove_tr=="yes" \
		and re.search(r'(\.[2-9]_[A-Z]{5})', seqfasta, flags=0)==None \
		and re.search(r'(\.[1-9]{2}_[A-Z]{5})', seqfasta, flags=0)==None \
		and re.search(r'(_T0[2-9]_[A-Z]{5})', seqfasta, flags=0)==None \
		and re.search(r'(_T[1-9]{2}_[A-Z]{5})', seqfasta, flags=0)==None \
		and re.search(r'(\.[2-9]-PROTEIN_[A-Z]{5})', seqfasta, flags=0)==None \
		and re.search(r'(\.[1-9]{2}-PROTEIN_[A-Z]{5})', seqfasta, flags=0)==None \
		and re.search(r'(\.[2-9]-PEP_[A-Z]{5})', seqfasta, flags=0)==None \
		and re.search(r'(\.[1-9]{2}-PEP_[A-Z]{5})', seqfasta, flags=0)==None \
		and re.search(r'(_00[2-9]_MUSAC)', seqfasta, flags=0)==None \
		and re.search(r'(_0[1-9]{2}_MUSAC)', seqfasta, flags=0)==None \
		and re.search(r'(_P0[2-9]_MAIZE)', seqfasta, flags=0)==None \
		and re.search(r'(_P[1-9]{2}_MAIZE)', seqfasta, flags=0)==None \
		and re.search(r'(g0[1-9]{2}_PHODC)', seqfasta, flags=0)==None \
		and re.search(r'(g00[2-9]_PHODC)', seqfasta, flags=0)==None : 
			output_file.write(seqfasta)
		elif remove_tr=="no":
			output_file.write(seqfasta)
		line = line.upper()
		line = line.replace("_ZEAMA","_MAIZE")
		line=line.replace('_ORYSA','_ORYSJ')
		line = line.replace("_SOLLY","_SOLLC")
		line = line.replace("_PHODA","_PHODC")
		line = line.replace('\"',"")
		spe_code=line.strip()[-5:]
		if re.search(r'(REF_)', line, flags=0) :
			spe_code = "REF_"+spe_code
# 		print spe_code
		locus_tag=line.strip()[1:-6]
		if spe_code in dic_of_dic.keys() and locus_tag in dic_of_dic[spe_code][spe_code].keys() :
			if "polypeptide_id" in dic_of_dic[spe_code][spe_code][locus_tag].keys():
				name = dic_of_dic[spe_code][spe_code][locus_tag]["polypeptide_id"]
				line = line.replace(locus_tag,name)
			elif "genfam_locus_tag" in dic_of_dic[spe_code][spe_code][locus_tag].keys():	
				name = locus_tag+".1"
				if name in dic_of_dic[spe_code][spe_code].keys() and "polypeptide_id" in dic_of_dic[spe_code][spe_code][locus_tag].keys():
					name = dic_of_dic[spe_code][spe_code][name]["polypeptide_id"]
					line = line.replace(locus_tag,name)
				else:
					name = dic_of_dic[spe_code][spe_code][locus_tag]["genfam_locus_tag"]
					name=name[:-6].replace("_GF","_PF")
					line = line.replace(locus_tag,name)
					
		# 	elif "genfam_locus_tag" in dic_of_dic[spe_code][spe_code][locus_tag].keys():
# 				name = dic_of_dic[spe_code][spe_code][locus_tag]["genfam_locus_tag"][:-6]
# 				line = line.replace(locus_tag,name)		
				
			
# 		if re.search(pat_vitvibis, line, flags=0):
# 			arathres=re.search(pat_vitvibis, line, flags=0)
# 			name = arathres.group(0).strip()
# 			newname=name+"_VITVI"
# 			line = line.replace(name,newname)
# 		if re.search(pat_vitvi, line, flags=0):
# 			arathres=re.search(pat_vitvi, line, flags=0)
# 			name = arathres.group(0).strip()
# 			newname=name+"_VITVI"
# 			line = line.replace(name,newname)
# 		elif re.search(pat_thecc, line, flags=0):
# 			arathres=re.search(pat_thecc, line, flags=0)
# 			name = arathres.group(0).strip()
# 			line = line.replace("t","g")
# 			newname=name+"_THECC"
# 			line = line.replace(name,newname)
# 		elif re.search(pat_theccbis, line, flags=0):
# 			arathres=re.search(pat_theccbis, line, flags=0)
# 			name = arathres.group(0).strip()
# 			newname=name+"_THECC"
# 			line = line.replace(name,newname)
# 		elif re.search(pat_orysj_MSU, line, flags=0): #si locus tag au format msu
# 			arathres=re.search(pat_orysj_MSU, line, flags=0)
# 			name = arathres.group(0).strip()
# 			if re.search("LOC_", name, flags=0):
# 				gene_name =name
# 			else:
# 				gene_name = "LOC_"+name
# 			if re.search("\.[0-9]+", name, flags=0):
# 				gene_name = gene_name
# 			else:
# 				gene_name = gene_name+".1"
# 			newname= dic_msu_to_rap[gene_name[:-1]+"1"]
# 			if newname!="None":
# 				newname=newname+gene_name[-2:]
# 			else:
# 				newname = gene_name
# 			if re.search("_ORYSJ", line, flags=0):
# 				newname=newname
# 			else :
# 				newname=newname+"_ORYSJ"
# 			line = line.replace(name,newname)
# 		elif re.search(pat_orysj_RAP, line, flags=0): #si locus tag au format RAP
# 			arathres=re.search(pat_orysj_RAP, line, flags=0)
# 			name = arathres.group(0).strip()
# 			if re.search("_ORYSJ", line, flags=0):
# 				newname=name
# 			else:
# 				if name[-1:]=="a":
# 					newname = name[0:-1]+".1"
# 				elif name[-1:]=="b":
# 					newname = name[0:-1]+".2"
# 				else:
# 					newname= name
# 				newname=newname+"_ORYSJ"
# 			line = line.replace(name,newname)
# 		elif re.search(pat_PHODC, line, flags=0):
# 			arathres=re.search(pat_PHODC, line, flags=0)
# 			name = arathres.group(0).strip()
# 			newname = name.replace("L","g")
# 			newname=newname+"_PHODC"
# 			line = line.replace(name,newname)
# 		elif re.search(pat_MALDO, line, flags=0):
# 			arathres=re.search(pat_MALDO, line, flags=0)
# 			name = arathres.group(0).strip()
# 			name = name.replace("C","P")
# 			newname=name+"_MALDO"
# 			line = line.replace(name,newname)
# 		elif re.search(pat_cucsa, line, flags=0):
# 			arathres=re.search(pat_cucsa, line, flags=0)
# 			name = arathres.group(0).strip()
# 			if re.search("_CUCSA", line, flags=0):
# 				newname=name
# 			else:
# 				newname=name+"_CUCSA"
# 			newname = newname.replace('CDS:','')
# 			line = line.replace(name,newname)
# 			line = line.replace(".1_","_")
# 			# print line
# 		elif re.search(pat_POPTR, line, flags=0):
# 			arathres=re.search(pat_POPTR, line, flags=0)
# 			name = arathres.group(0)
# 			newname=name[:-1]+"_POPTR"
# 			line = line.replace(name,newname)
# 		elif re.search(pat_MUSBA, line, flags=0):
# 			arathres=re.search(pat_MUSBA, line, flags=0)
# 			name = arathres.group(0)
# 			newname=name[:-1]+"_MUSBA"
# 			newname = newname.replace('G','T')
# 			newname = newname.replace('P','T')
# 			line = line.replace(name,newname)
# 		elif re.search(pat_CAJCA, line, flags=0):
# 			arathres=re.search(pat_CAJCA, line, flags=0)
# 			name = arathres.group(0)
# 			newname=name[:-1]+"_CAJCA"
# 			line = line.replace(name,newname)
# 		elif re.search(pat_horvu, line, flags=0):
# 			arathres=re.search(pat_horvu, line, flags=0)
# 			name = arathres.group(0)
# 			newname=name[:-1]+"_HORVU"
# 			line = line.replace(name,newname)
# 		elif re.search(pat_citsi, line, flags=0):
# 			arathres=re.search(pat_citsi, line, flags=0)
# 			name = arathres.group(0)
# 			newname=name[:-1]+"_CITSI"
# 			line = line.replace(name,newname)
# 		elif re.search(pat_sorbi, line, flags=0):
# 			arathres=re.search(pat_sorbi, line, flags=0)
# 			name = arathres.group(0)
# 			newname=name[:-1]+"_SORBI"
# 			line = line.replace(name,newname)
# 		elif re.search(pat_soltu, line, flags=0):
# 			arathres=re.search(pat_soltu, line, flags=0)
# 			name = arathres.group(0)
# 			newname=name[:-1]+"_SOLTU"
# 			newname = newname.replace('T','G')
# 			newname = newname.replace('P','G')
# 			line = line.replace(name,newname)
# 		elif re.search(pat_sollc, line, flags=0):
# 			arathres=re.search(pat_sollc, line, flags=0)
# 			name = arathres.group(0)
# 			newname=name[:-1]+"_SOLLC"
# 			line = line.replace(name,newname)
# 		elif re.search(pat_ricco, line, flags=0):
# 			arathres=re.search(pat_ricco, line, flags=0)
# 			name = arathres.group(0)
# 			newname=name[:-1]+"_RICCO"
# 			line = line.replace(name,newname)
# 		elif re.search(pat_ARATH_transc, line, flags=0):
# 			arathres=re.search(pat_ARATH_transc, line, flags=0)
# 			name = arathres.group(0)
# 			newname=name[:-1]+"_ARATH"
# 			line = line.replace(name,newname)
# 		elif re.search(pat_ARATH, line, flags=0):
# 			arathres=re.search(pat_ARATH, line, flags=0)
# 			name = arathres.group(0)[:-1]
# 			newname=name+".1_ARATH"
# 			line = line.replace(name,newname)
# 		elif re.search(pat_gosra, line, flags=0):
# 			arathres=re.search(pat_gosra, line, flags=0)
# 			name = arathres.group(0)
# 			newname=name[:-1]+"_GOSRA"
# 			line = line.replace(name,newname)
# 		elif re.search(pat_lotja, line, flags=0):
# 			arathres=re.search(pat_gosra, line, flags=0)
# 			name = arathres.group(0)
# 			newname=name[:-1]+"_LOTJA"
# 			line = line.replace(name,newname)
# 		elif re.search(pat_musa, line, flags=0):
# 			arathres=re.search(pat_musa, line, flags=0)
# 			name = arathres.group(0)[:-1]
# 			newname = name.replace('T','G')
# 			newname = newname.replace('P','G')
# 			newname=name+"_MUSAC"
# 			line = line.replace(name,newname)
# 		elif re.search(pat_medtr, line, flags=0):
# 			arathres=re.search(pat_medtr, line, flags=0)
# 			name = arathres.group(0)
# 			newname=name[:-1]+".1_MEDTR"
# 			line = line.replace(name,newname)
# 		elif re.search(pat_medtr_transcr, line, flags=0):
# 			arathres=re.search(pat_medtr_transcr, line, flags=0)
# 			name = arathres.group(0)
# 			if re.search("_MEDTR", line, flags=0):
# 				newname=name
# 			else:
# 				newname=name+"_MEDTR"
# 			line = line.replace(name,newname)
# 		elif re.search(pat_manes, line, flags=0):
# 			manesres=re.search(pat_manes, line, flags=0)
# 			name = manesres.group(0)
# 			newname=name[:-1]+"_MANES"
# 			line = line.replace(name,newname)
# 		elif re.search(pat_maize, line, flags=0):
# 			bradires=re.search(pat_maize, line, flags=0)
# 			name = bradires.group(0)[:-1]
# 			newname=name.replace(".v6a","_MAIZE")
# 			newname=newname.replace("P","G")
# 			newname=newname.replace("T","G")
# 			newname=name+"_MAIZE"
# 			line = line.replace(name,newname)
# 		elif re.search(pat_glyma, line, flags=0):
# 			bradires=re.search(pat_glyma, line, flags=0)
# 			name = bradires.group(0)
# 			newname=name[:-1]+".1_GLYMA"
# 			line = line.replace(name,newname)
# 		elif re.search(pat_glyma_transc, line, flags=0):
# 			bradires=re.search(pat_glyma_transc, line, flags=0)
# 			name = bradires.group(0)
# 			newname=name[:-1]+"_GLYMA"
# 			line = line.replace(name,newname)
# 		elif re.search(pat_bradi_transc, line, flags=0):
# 			bradires=re.search(pat_bradi_transc, line, flags=0)
# 			name = bradires.group(0)
# 			newname=name[:-1]+"_BRADI"
# 			line = line.replace(name,newname)
# 		elif re.search(pat_bradi, line, flags=0):
# 			bradires=re.search(pat_bradi, line, flags=0)
# 			name = bradires.group(0)
# 			newname=name[:-1]+".1_BRADI"
# 			line = line.replace(name,newname)
		seqfasta=line
	else:
		seqfasta=seqfasta+line

if remove_tr=="yes" \
and re.search(r'(\.[2-9]_[A-Z]{5})', seqfasta, flags=0)==None \
and re.search(r'(\.[1-9]{2}_[A-Z]{5})', seqfasta, flags=0)==None \
and re.search(r'(_T0[2-9]_[A-Z]{5})', seqfasta, flags=0)==None \
and re.search(r'(_T[1-9]{2}_[A-Z]{5})', seqfasta, flags=0)==None \
and re.search(r'(\.[2-9]-PROTEIN_ARATH)', seqfasta, flags=0)==None \
and re.search(r'(\.[1-9]{2}-PROTEIN_ARATH)', seqfasta, flags=0)==None \
and re.search(r'(\.[2-9]-PEP_[A-Z]{5})', seqfasta, flags=0)==None \
and re.search(r'(\.[1-9]{2}-PEP_[A-Z]{5})', seqfasta, flags=0)==None \
and re.search(r'(_00[2-9]_MUSAC)', seqfasta, flags=0)==None \
and re.search(r'(_0[1-9]{2}_MUSAC)', seqfasta, flags=0)==None \
and re.search(r'(_P0[2-9]_MAIZE)', seqfasta, flags=0)==None \
and re.search(r'(_P[1-9]{2}_MAIZE)', seqfasta, flags=0)==None \
and re.search(r'(g0[1-9]{2}_PHODC)', seqfasta, flags=0)==None \
and re.search(r'(g00[2-9]_PHODC)', seqfasta, flags=0)==None :
	output_file.write(seqfasta)
elif remove_tr=="no":
	output_file.write(seqfasta)	

			
fasta_file.close()
output_file.close()

