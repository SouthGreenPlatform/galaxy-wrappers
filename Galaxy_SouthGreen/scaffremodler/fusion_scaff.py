
#
#  Copyright 2014 CIRAD
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, see <http://www.gnu.org/licenses/> or
#  write to the Free Software Foundation, Inc.,
#  51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#
#

import optparse, os, shutil, subprocess, sys, tempfile, fileinput, ConfigParser, operator, time, random


from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def rev_seq(seq):
	#function that reverse and complement a sequence
	my_dna = Seq(seq, generic_dna)
	return str(my_dna.reverse_complement())

def charge_data(TABLE):
	#function that record scaffold target and destination
	dic_dest = {}#hash table : key = destination name, word = [destination orientation, last position of destination, hash_tableX, liste of destination coordinates]
	#hash_tableX : key = last position of destination, word = [(0) target_name, (1) target_beginning, (2) target_end, (3) target_orientation, (4) destination_name, (5) destination_beginning, (6) destination_end, (7) destination_orientation, (8) Fusion_type]
	dic_cible = set()#record scaffold target
	file = open(TABLE)
	for line in file:
		data = line.split()
		if data != []:
			if data[0] in dic_dest:
				sys.exit('Problem in the table file target is also destination: '+data[0])
			elif data[4] in dic_cible:
				sys.exit('Problem in the table file destination is also a target: '+data[4])
			dic_cible.add(data[0])
			if data[4] in dic_dest:
				if dic_dest[data[4]][0] != data[7]:
					sys.exit('Problem in the table file destination orientation is changing : '+data[4])
				else:
					dic_dest[data[4]][2][int(data[6])] = [data[0], str(int(data[1])-1), str(int(data[2])-1), data[3], data[4], str(int(data[5])-1), str(int(data[6])-1), data[7], data[8]]
					dic_dest[data[4]][1].append((int(data[6])-1))
					if [str(int(data[5])-1), str(int(data[6])-1)] in dic_dest[data[4]][3]:
						sys.exit('Problem in the table file. There is more than one target for same position in destination : '+data[4]+' '+data[5]+' '+data[6])
					dic_dest[data[4]][3].append([str(int(data[5])-1), str(int(data[6])-1)])
			else:
				dic_dest[data[4]] = [data[7], [(int(data[6])-1)],{},[]]
				dic_dest[data[4]][2][(int(data[6])-1)] = [data[0], str(int(data[1])-1), str(int(data[2])-1), data[3], data[4], str(int(data[5])-1), str(int(data[6])-1), data[7], data[8]]
				dic_dest[data[4]][3].append([str(int(data[5])-1), str(int(data[6])-1)])
	return dic_dest

def verif_borne(DEBUT, FIN, SEQ, NOM):
	# os.system("echo '"+NOM+" "+str(DEBUT)+" "+str(FIN)+"'")
	if DEBUT > FIN:
		sys.exit('Problem (1) in the table file wrong target coordinates : '+NOM+' '+str(DEBUT)+' '+str(FIN))
	if FIN + 1 > len(SEQ):
		sys.exit('Problem (2) in the table file wrong target coordinates : '+NOM+' '+str(FIN))
	if SEQ[DEBUT] == 'n' or SEQ[DEBUT] == 'N':
		sys.exit('Problem (3) in the table file wrong target coordinates : '+NOM+' '+str(DEBUT))
	if SEQ[FIN] == 'n' or SEQ[FIN] == 'N':
		sys.exit('Problem (4) in the table file wrong target coordinates : '+NOM+' '+str(FIN)+' '+SEQ[FIN-1])
	if len(SEQ) != (FIN - DEBUT + 1):
		if DEBUT != 0:
			if SEQ[DEBUT-1] != 'n' and SEQ[DEBUT-1] != 'N':
				sys.exit('Problem (5) in the table file wrong target coordinates : '+NOM+' '+str(DEBUT))
		if FIN + 1 != len(SEQ):
			if SEQ[FIN+1] != 'n' and SEQ[FIN+1] != 'N':
				sys.exit('Problem (6) in the table file wrong target coordinates : '+NOM+' '+str(FIN))

def verif_no_nuc(TEST, START, END, ID):
	if START != END:
		TEST2 = TEST[START+1:END-1]
		if 'a' in TEST2 or 'A' in TEST2 or 't' in TEST2 or 'T' in TEST2 or 'g' in TEST2 or 'G' in TEST2 or 'c' in TEST2 or 'C' in TEST2:
			sys.exit('There are nucleotides in the sequence that is considered as unknown in '+ID+' '+str(START)+' '+str(END))
		elif 'X' in TEST2:
			sys.exit('There are X in the sequence that are considered as unknown in '+ID+' '+str(START)+' '+str(END))
	

def reconstruct(DIC, LISTE, OR, DEST, DICTIO, OUT):
	#0)creation de la liste des cibles et destination 1)recuperer la destination 2)recuperer les cibles
	#os.system("echo '***********loading sequence**********'")
	dico = {}
	dico[DEST] = str(DICTIO[DEST].seq)
	dico_or = {}
	dico_or[DEST] = OR
	LISTE_TAR = {}
	for n in DIC:
		if DIC[n][0] in LISTE_TAR:
			LISTE_TAR[DIC[n][0]].append([int(DIC[n][1]),int(DIC[n][2])])
		else:
			LISTE_TAR[DIC[n][0]] = [[int(DIC[n][1]),int(DIC[n][2])]]
		dico[DIC[n][0]] = str(DICTIO[DIC[n][0]].seq)
		if DIC[n][0] in dico_or:#verify orientation of target if several blocs
			if dico_or[DIC[n][0]] != DIC[n][3]:
				sys.exit('Problem in the table file target orientation is changing : '+DIC[n][0])
			dico_or[DIC[n][0]] = DIC[n][3]
		else:
			dico_or[DIC[n][0]] = DIC[n][3]
		#Additionale verification for destination
		if OR == 'REV':
			verif_no_nuc(rev_seq(dico[DEST]), int(DIC[n][5]), int(DIC[n][6]), DEST)
		else:
			verif_no_nuc(dico[DEST], int(DIC[n][5]), int(DIC[n][6]), DEST)
	#Additionale verification for target
	for n in LISTE_TAR:
		liste_tar = sorted(LISTE_TAR[n], key=operator.itemgetter(0))
		if liste_tar[0][0] != 0:
			sys.exit('Problem in the table file it is not the full size target (the begining is missing) : '+n+' '+str(liste_tar[0][0]))
		while len(liste_tar) > 1:
			verif_no_nuc(dico[n], liste_tar[0][1], liste_tar[1][0], n)
			del liste_tar[0]
		if len(dico[n]) != liste_tar[0][1]+1:
			sys.exit('Problem in the table file it is not the full size target (the end is missing) : '+n+' '+str(liste_tar[0][1]))
	#os.system("echo 'sequence loaded'")
	#for n in dico:
		#os.system("echo '"+str(n)+"\t"+str(len(dico[n]))+"'")
	#3)les orienter en fonction de l'info donnee
	#os.system("echo 'orienting sequence'")
	for n in dico_or:
		if dico_or[n] == 'REV':
			dico[n] = rev_seq(dico[n])
		#os.system("echo '"+n+" "+dico_or[n]+"'")
	#os.system("echo 'sequence orientated'")
	#4)reconstruire de super_scaffold
	#os.system("echo 'recontructing sequence'")
	fin = len(dico[DEST])
	sequence = ''
	liste_scaff = []
	liste_pos = []
	# print LISTE
	# print DIC
	for n in LISTE:
		if int(DIC[n][5]) == int(DIC[n][6]) and int(DIC[n][6]) != 0:#there is a sequence at the end of the destination
			if fin == int(DIC[n][5]) + 1:
				liste_pos = calcul_new_pos(liste_pos, (int(DIC[n][2])+1 - int(DIC[n][1])), 20)
				# print liste_pos
				liste_scaff.insert(0, str(DIC[n][0]))
				verif_borne(int(DIC[n][1]), int(DIC[n][2]), dico[DIC[n][0]], DIC[n][0])
				if DIC[n][8] == 'REV':
					sequence = rev_seq(dico[DIC[n][0]][int(DIC[n][1]):int(DIC[n][2])+1]) + sequence
				else:
					sequence = dico[DIC[n][0]][int(DIC[n][1]):int(DIC[n][2])+1] + sequence
				sequence = 'nnnnnnnnnnnnnnnnnnnn' + sequence
			else:
				sys.exit('Problem in the table file destination expected at the end but does not correspond: '+DIC[n][0])
		else:
			liste_scaff.insert(0, DEST)
			verif_borne(int(DIC[n][6]), fin-1, dico[DEST], DEST)
			sequence = dico[DEST][int(DIC[n][6]):fin] + sequence
			if int(DIC[n][5]) == 0 and int(DIC[n][6]) == 0:#there is sequence at the beginning of the destination
				sequence = 'nnnnnnnnnnnnnnnnnnnn' + sequence
				liste_pos = calcul_new_pos(liste_pos, (fin - int(DIC[n][6])), 20)
			else:#there is NO sequence at the beginning of the destination
				if int((((int(DIC[n][6]) - int(DIC[n][5])) - 1)-(int(DIC[n][2])+1 - int(DIC[n][1])))/2) > 0:
					liste_pos = calcul_new_pos(liste_pos, (fin - int(DIC[n][6])), int((((int(DIC[n][6]) - int(DIC[n][5])) - 1)-(int(DIC[n][2])+1 - int(DIC[n][1])))/2))
					j = 0
					while j < int((((int(DIC[n][6]) - int(DIC[n][5])) - 1)-(int(DIC[n][2])+1 - int(DIC[n][1])))/2):
						sequence = 'n' + sequence
						j = j + 1
					j = 0
				else:
					sequence = 'nnnnnnnnnnnnnnnnnnnn' + sequence
					liste_pos = calcul_new_pos(liste_pos, (fin - int(DIC[n][6])), 20)
			# print liste_pos
			fin = int(DIC[n][5])+1
			liste_scaff.insert(0, str(DIC[n][0]))
			verif_borne(int(DIC[n][1]), int(DIC[n][2]), dico[DIC[n][0]], DIC[n][0])
			if DIC[n][8] == 'REV':
				sequence = rev_seq(dico[DIC[n][0]][int(DIC[n][1]):int(DIC[n][2])+1]) + sequence
			else:
				sequence = dico[DIC[n][0]][int(DIC[n][1]):int(DIC[n][2])+1] + sequence
			if int(DIC[n][5]) != 0 and int(DIC[n][6]) != 0:#there is NO sequence at the beginning of the destination
				if int((((int(DIC[n][6]) - int(DIC[n][5])) - 1)-(int(DIC[n][2])+1 - int(DIC[n][1])))/2) > 0:
					liste_pos = calcul_new_pos(liste_pos, (int(DIC[n][2])+1 - int(DIC[n][1])), int((((int(DIC[n][6]) - int(DIC[n][5])) - 1)-(int(DIC[n][2])+1 - int(DIC[n][1])))/2))
					j = 0
					while j < int((((int(DIC[n][6]) - int(DIC[n][5])) - 1)-(int(DIC[n][2])+1 - int(DIC[n][1])))/2):
						sequence = 'n' + sequence
						j = j + 1
					j = 0
				else:
					sequence = 'nnnnnnnnnnnnnnnnnnnn' + sequence
					liste_pos = calcul_new_pos(liste_pos, (int(DIC[n][2])+1 - int(DIC[n][1])), 20)
			else:#there is sequence at the beginning of the destination
				liste_pos = calcul_new_pos(liste_pos, (int(DIC[n][2])+1 - int(DIC[n][1])), 0)
			# print liste_pos
	if int(DIC[n][5]) != 0 and int(DIC[n][6]) != 0:#there is NO sequence at the beginning of the destination
		sequence = dico[DEST][0:fin] + sequence
		liste_pos = calcul_new_pos(liste_pos, fin, 0)
		liste_scaff.insert(0, DEST)
		verif_borne(0, fin-1, dico[DEST], DEST)
		# print liste_pos
	#print '******************************', liste_pos
	#os.system("echo 'lenght seq final :"+str(len(sequence))+"'")
	#os.system("echo '**********sequence reconstructed**********'")
	#5)retourner le scaffold dans la bonne orientation
	mot = ''
	if OR == 'REV':
		i = 0
		for n in liste_scaff:
			if i == 0:
				mot = n+'\t'+str((abs(liste_pos[i+1]-(len(sequence)-1)))+1)+'\t'+str((abs(liste_pos[i]-(len(sequence)-1)))+1)
			else:
				mot = n+'\t'+str((abs(liste_pos[i+1]-(len(sequence)-1)))+1)+'\t'+str((abs(liste_pos[i]-(len(sequence)-1)))+1)+'\t'+mot
			i = i + 2
		outfile = open(OUT,'a')
		outfile.write(DEST+'\t'+mot+'\n')
		outfile.close()
		# print len(sequence)
		#os.system("echo 'reversing sequence'")
		return rev_seq(sequence)
	else:
		i = 0
		for n in liste_scaff:
			if i == 0:
				mot = n+'\t'+str(liste_pos[i]+1)+'\t'+str(liste_pos[i+1]+1)
			else:
				mot = mot+'\t'+n+'\t'+str(liste_pos[i]+1)+'\t'+str(liste_pos[i+1]+1)
			i = i + 2
		outfile = open(OUT,'a')
		outfile.write(DEST+'\t'+mot+'\n')
		outfile.close()
		# print len(sequence)
		return sequence

def calcul_new_pos(l1, taille_seq, ADD):
	if len(l1) == 0:
		l1 = [ADD, ADD + taille_seq-1]
		return l1
	else:
		LISTE_A_RET = []
		for n in l1:
			LISTE_A_RET.append(n+ADD+taille_seq)
		LISTE_A_RET.insert(0,ADD+taille_seq-1)
		LISTE_A_RET.insert(0,ADD)
		return LISTE_A_RET

def __main__():
	#Parse Command Line
	parser = optparse.OptionParser(usage="python %prog [options]\n\nProgram designed by Guillaume MARTIN : guillaume.martin@cirad.fr\n\n"
	"This script merge scaffold based on tabulated file")
	# Wrapper options.
	parser.add_option( '', '--table', dest='table', default='not_filled', help='The table file of scaffold to merge')
	parser.add_option( '', '--fasta', dest='fasta', default='not_filled', help='The multi-fasta scaffold file')
	parser.add_option( '', '--out', dest='out', default='fusion.fasta', help='The multi-fasta output file name')
	parser.add_option( '', '--out_verif', dest='out_verif', default='fusion2verif.txt', help='The output file to give to verif_fusion.py')
	(options, args) = parser.parse_args()
	

	dico = charge_data(options.table)


	record_dict = SeqIO.index(options.fasta, "fasta")
	
	outfile = open(options.out_verif,'w')
	outfile.close()
	
	outfile = open(options.out,'w')
	for j in dico:
		liste = sorted(dico[j][1], reverse=True)
		SEQ = SeqRecord(Seq(reconstruct(dico[j][2], liste, dico[j][0], j, record_dict, options.out_verif), generic_dna), id = j, description='')
		SeqIO.write(SEQ, outfile, "fasta")

	file = open(options.table)
	dico = set()
	for line in file:
		data = line.split()
		if data != []:
			dico.add(data[0])
			dico.add(data[4])

	for n in record_dict:
		if not(n in dico):
			SeqIO.write(record_dict[n], outfile, "fasta")
	outfile.close()

if __name__ == "__main__": __main__()