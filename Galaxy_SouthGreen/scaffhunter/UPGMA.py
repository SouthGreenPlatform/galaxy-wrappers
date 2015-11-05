#!/usr/local/bioinfo/python/2.7.9/bin/python
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

def cherche_max_couple(DICO, INDEX):
	maxi  = 0
	i = 0
	taille = len(INDEX)
	while i < taille:
		for n in DICO:
			if n != INDEX[i]:
				if DICO[n][i] > maxi:
					maxi = DICO[n][i]
					couple = [n,INDEX[i],((1 - maxi/100)/2)]
		i += 1
	return couple

def cherche_min_couple(DICO, INDEX):
	mini  = 10000000000000
	i = 0
	taille = len(INDEX)
	while i < taille:
		for n in DICO:
			if n != INDEX[i]:
				if DICO[n][i] < mini:
					mini = DICO[n][i]
					couple = [n,INDEX[i],(mini/2)]
		i += 1
	return couple


def pseudoUP(MAT, IND, INFO, TYPE, MATRIX_MARK, LISTE_ORDRE, DICO_MARK, DICO_ORIENT, DICO_INDEX, NB_MARK):
	index_mat = list(IND)
	dico = {}
	dico_couple = {}
	for n in MAT:
		dico[n] = list(MAT[n])
	i = 0
	while len(dico) > 1:
		i += 1
		#On cherche le meilleur couple
		if TYPE == 'IDENT':
			couple = cherche_max_couple(dico, index_mat)
		elif TYPE == 'DIF':
			couple = cherche_min_couple(dico, index_mat)
		else:
			sys.exit('Problem in data type')
		dico_couple['g'+str(i)] = couple
		#####On recalcule la table
		in_line = []
		j = 0
		n1 = INFO[couple[0]]
		n2 = INFO[couple[1]]
		while j < len(index_mat):
			val1 = dico[couple[0]][j]
			val2 = dico[couple[1]][j]
			in_line.append(((val1*n1)+(val2*n2))/(n1+n2))
			j += 1
		in_line.append(100)
		##On enregistre les nouvelles infos
		INFO['g'+str(i)] = n1+n2
		dico['g'+str(i)] = in_line
		j = 0
		for n in index_mat:
			dico[n].append(in_line[j])
			j += 1
		index_mat.append('g'+str(i))##tres important que ca arrive apres la boucle precedente sinon la derniere case du tableau est duplique
		##on retire ce qui est en trop
		if index_mat.index(couple[0]) < index_mat.index(couple[1]):
			c1 = index_mat.index(couple[0])
			c2 = index_mat.index(couple[1])
		else:
			c1 = index_mat.index(couple[1])
			c2 = index_mat.index(couple[0])
		#On retire les colonnes
		for n in dico:
			del dico[n][c2]
			del dico[n][c1]
		#On retire les lignes
		del dico[couple[0]]
		del dico[couple[1]]
		#On retire aussi les indexes et les info de nombre de marqueurs
		del index_mat[c2]
		del index_mat[c1]
		del INFO[couple[0]]
		del INFO[couple[1]]
		##Une petite verification
		if couple[0] in INFO or couple[1] in INFO or couple[0] in dico or couple[1] in dico or len(dico) != len(INFO) or len(INFO) != len(index_mat):
			sys.exit('There is a bug')
		taille = len(dico)
		for n in dico:
			if len(dico[n]) != taille:
				print len(dico[n]), taille
				sys.exit('There is a bug')
	#Creation of the tree structure
	return create_nexus(dico_couple, MAT, IND, MATRIX_MARK, LISTE_ORDRE, DICO_MARK, DICO_ORIENT, DICO_INDEX, NB_MARK, TYPE)

def create_nexus(DICO, MAT, IND, MATRIX_MARK, LISTE_ORDRE, DICO_MARK, DICO_ORIENT, DICO_INDEX, NB_MARK, TYPE):
	dico = {}
	for n in DICO:
		dico[n] = list(DICO[n])
	taille = len(dico)
	i = 0
	dico_mot = {}
	dico_arbre = {}
	while i < taille:
		i += 1
		if dico['g'+str(i)][0] in dico_arbre and dico['g'+str(i)][1] in dico_arbre:
			dico_arbre['g'+str(i)] = '('+dico_arbre[dico['g'+str(i)][0]]+':'+str(dico['g'+str(i)][2])+','+dico_arbre[dico['g'+str(i)][1]]+':'+str(dico['g'+str(i)][2])+')'
			dico_mot['g'+str(i)] = cherche_best_pos(dico_mot[dico['g'+str(i)][0]], dico_mot[dico['g'+str(i)][1]], MAT, IND, MATRIX_MARK, LISTE_ORDRE, DICO_MARK, DICO_ORIENT, DICO_INDEX, NB_MARK, TYPE)
		elif dico['g'+str(i)][0] in dico_arbre:
			dico_arbre['g'+str(i)] = '('+dico_arbre[dico['g'+str(i)][0]]+':'+str(dico['g'+str(i)][2])+','+dico['g'+str(i)][1]+':'+str(dico['g'+str(i)][2])+')'
			dico_mot['g'+str(i)] = cherche_best_pos(dico_mot[dico['g'+str(i)][0]], [dico['g'+str(i)][1]], MAT, IND, MATRIX_MARK, LISTE_ORDRE, DICO_MARK, DICO_ORIENT, DICO_INDEX, NB_MARK, TYPE)
		elif dico['g'+str(i)][1] in dico_arbre:
			dico_arbre['g'+str(i)] = '('+dico['g'+str(i)][0]+':'+str(dico['g'+str(i)][2])+','+dico_arbre[dico['g'+str(i)][1]]+':'+str(dico['g'+str(i)][2])+')'
			dico_mot['g'+str(i)] = cherche_best_pos(dico_mot[dico['g'+str(i)][1]], [dico['g'+str(i)][0]], MAT, IND, MATRIX_MARK, LISTE_ORDRE, DICO_MARK, DICO_ORIENT, DICO_INDEX, NB_MARK, TYPE)
		else:
			dico_arbre['g'+str(i)] = '('+dico['g'+str(i)][0]+':'+str(dico['g'+str(i)][2])+','+dico['g'+str(i)][1]+':'+str(dico['g'+str(i)][2])+')'
			dico_mot['g'+str(i)] = cherche_best_pos([dico['g'+str(i)][0]], [dico['g'+str(i)][1]], MAT, IND, MATRIX_MARK, LISTE_ORDRE, DICO_MARK, DICO_ORIENT, DICO_INDEX, NB_MARK, TYPE)
	return [dico_arbre['g'+str(i)], dico_mot['g'+str(i)]]

def cherche_best_pos(G1, G2, MAT, IND, MATRIX_MARK, LISTE_ORDRE, DICO_MARK, DICO_ORIENT, DICO_INDEX, NB_MARK, TYPE):
	ordre1 = G1+G2
	ordre2 = G2+G1
	reverse2 = list(G2)
	reverse2.reverse()
	ordre3 = G1+reverse2
	ordre4 = reverse2+G1
	dico_reorient = {}
	for n in DICO_ORIENT:
		if n in G2:
			if DICO_ORIENT[n] == 'FWD':
				dico_reorient[n] = 'REV'
			elif DICO_ORIENT[n] == 'REV':
				dico_reorient[n] = 'FWD'
			else:
				sys.exit('bug in cherche_best_pos1')
		else:
			dico_reorient[n] = str(DICO_ORIENT[n])
	if TYPE == 'IDENT':
		score1 = calcul_score_id(MATRIX_MARK, ordre1, DICO_MARK, DICO_ORIENT, DICO_INDEX, NB_MARK)
		score2 = calcul_score_id(MATRIX_MARK, ordre2, DICO_MARK, DICO_ORIENT, DICO_INDEX, NB_MARK)
		score3 = calcul_score_id(MATRIX_MARK, ordre3, DICO_MARK, dico_reorient, DICO_INDEX, NB_MARK)
		score4 = calcul_score_id(MATRIX_MARK, ordre4, DICO_MARK, dico_reorient, DICO_INDEX, NB_MARK)
	elif TYPE == 'DIF':
		score1 = calcul_score_dif(MATRIX_MARK, ordre1, DICO_MARK, DICO_ORIENT, DICO_INDEX, NB_MARK)
		score2 = calcul_score_dif(MATRIX_MARK, ordre2, DICO_MARK, DICO_ORIENT, DICO_INDEX, NB_MARK)
		score3 = calcul_score_dif(MATRIX_MARK, ordre3, DICO_MARK, dico_reorient, DICO_INDEX, NB_MARK)
		score4 = calcul_score_dif(MATRIX_MARK, ordre4, DICO_MARK, dico_reorient, DICO_INDEX, NB_MARK)
	else:
		sys.exit('bug in cherche_best_pos2')
	liste = [score1, score2, score3, score4]
	INDEX = liste.index(max(liste))
	if INDEX == 0:
		return list(ordre1)
	elif INDEX == 1:
		return list(ordre2)
	elif INDEX == 2:
		for n in dico_reorient:
			DICO_ORIENT[n] = str(dico_reorient[n])
		return list(ordre3)
	elif INDEX == 3:
		for n in dico_reorient:
			DICO_ORIENT[n] = str(dico_reorient[n])
		return list(ordre4)
	else:
		sys.exit('bug in cherche_best_pos3')

def calcul_score_id(MAT, ORDRE, MARK, ORIENT, DIC_INDEX, NB):
	#on cree l'ordre des markers
	scaff_ordre = list(ORDRE)
	ordre = []
	for n in scaff_ordre:
		liste = list(MARK[n])
		if ORIENT[n] == 'REV':
			liste.reverse()
		ordre = ordre + liste
	#on enregistre les positions relative des markers dans cet ordre
	liste_ordre = {}
	i = 0
	for n in ordre:
		liste_ordre[n] = i
		i += 1
	#on calcule le score
	score = 0
	while len(scaff_ordre) > 1:
		scaff_fait = scaff_ordre[0]
		del scaff_ordre[0]
		# print scaff_fait, MARK[scaff_fait]
		for n in MARK[scaff_fait]:
			mark = DIC_INDEX[n]
			for j in scaff_ordre:
				for k in MARK[j]:
					val = MAT[mark][k]
					if val != 999999999:
						score = score + ((1-(abs(liste_ordre[n]-liste_ordre[k])/NB))*val)
	return score
	
def calcul_score_dif(MAT, ORDRE, MARK, ORIENT, DIC_INDEX, NB):
	#on cree l'ordre des markers
	scaff_ordre = list(ORDRE)
	ordre = []
	for n in scaff_ordre:
		liste = list(MARK[n])
		if ORIENT[n] == 'REV':
			liste.reverse()
		ordre = ordre + liste
	#on enregistre les positions relative des markers dans cet ordre
	liste_ordre = {}
	i = 0
	for n in ordre:
		liste_ordre[n] = i
		i += 1
	#on calcule le score
	score = 0
	while len(scaff_ordre) > 1:
		scaff_fait = scaff_ordre[0]
		del scaff_ordre[0]
		# print scaff_fait, MARK[scaff_fait]
		for n in MARK[scaff_fait]:
			mark = DIC_INDEX[n]
			for j in scaff_ordre:
				for k in MARK[j]:
					val = MAT[mark][k]
					if val != 999999999:
						score = score + ((abs(liste_ordre[n]-liste_ordre[k])/NB)*val)
	return score

def record_scaff(FILE, HEAD):
	file = open(FILE)
	dico = {}
	for line in file:
		data = line.split()
		if data:
			if '#' in data[1]:
				sys.exit('Please change scaffold names otherwise there will be somme problemes. Removing "#" will be enough')
			if data[1] in dico:
				dico[data[1]][0].append(HEAD.index(data[0]))
				dico[data[1]][1].append(data[0])
			else:
				dico[data[1]] = [[],[]]
				dico[data[1]][0].append(HEAD.index(data[0]))
				dico[data[1]][1].append(data[0])
	file.close()
	return dico

def appar(MAT, INDEX, SCAFF, GROUP):
	dic = set()
	dico_couple = {}
	# t0 = time.clock()
	for n in SCAFF:
		for j in SCAFF:
			if not(j) in dic:
				liste_val = 0
				nb_val = 0
				for l in SCAFF[n][1]:
					for m in SCAFF[j][0]:
						val = MAT[l][m]
						if val != 999999999:
							liste_val = liste_val + MAT[l][m]
							nb_val += 1
				if nb_val:
					dico_couple[n+"#"+j] = str(liste_val/float(nb_val))
				else:
					dico_couple[n+"#"+j] = str(0)
		dic.add(n)
		# os.system('echo "'+n+' '+str(time.clock() - t0)+'"')
		# t0 = time.clock()
	return dico_couple

def mat_mark2mat_scaff(SCAFF_FILE, MAT, OUT):
	#on enregistre dans un dico les scaffolds et l'ordre des marqueurs dans le scaffold
	file = open(MAT)
	index_mat = file.readline().split()[1:]
	file.close()
	scaffold = record_scaff(SCAFF_FILE, index_mat)
	
	#Recording pairwise information
	#on charge la matrice dans un dico
	matrix = {}
	file = open(MAT)
	file.readline()
	for line in file:
		data = line.split()
		if data:
			if data[0] in matrix:
				sys.exit("Warning there is redundancy in markers names, verify your data")
			else:
				matrix[data[0]] = map(float, data[1:])
	file.close()
	
	#Calculating mean pairwise between scaffolds
	APAR = appar(matrix, index_mat, scaffold, [])
	
	#Creating the table:
	dico_table = {}
	index_table = []
	for n in scaffold:
		index_table.append(n)
	for n in index_table:
		dico_table[n] = []
		for j in index_table:
			dico_table[n].append('-')
	#Filling the table
	for n in APAR:
		scaff = n.split('#')
		dico_table[scaff[0]][index_table.index(scaff[1])] = APAR[n]
		dico_table[scaff[1]][index_table.index(scaff[0])] = APAR[n]
	
	#Printing results
	OUT.write('ID\tnb_mark\t'+'\t'.join(index_table)+'\n')
	for n in index_table:
		OUT.write(n+'\t'+str(len(scaffold[n][0]))+'\t'+'\t'.join(dico_table[n])+'\n')
	OUT.flush()

def __main__():
	#Parse Command Line
	parser = optparse.OptionParser(usage="python %prog [options]\n\nProgram designed by Guillaume MARTIN : guillaume.martin@cirad.fr\n\n"
	"This program takes in input a matrix containing maker linkage or divergence and a tabulated file locating these markers on scaffolds "
	"(or sequence) and try to group and order them based on an UPGMA like approach.\n"
	"This program work as followed:\n" 
	"\t1 - an mean linkage/divergence is calculated between scaffolds\n"
	"\t2 - scaffolds are then grouped using an UPGMA like approach\n"
	"\t3 - scaffolds are orientated and positioned (at the beginning or the end of a precedent group) trying to optimize a score calculated for each position and orientation.\n"
	"This program output a file containing ordered scaffolds in a table format and a table file containing ordered markers relative to the scaffold order.\n")
	# Wrapper options.
	parser.add_option( '', '--scaff', dest='scaff', default=None, help='A table file containing in column 1: markers name, column 2: scaffold name. Columns should be ordered first on scaffold name and second on marker position on scaffold')
	parser.add_option( '', '--mat', dest='mat', default=None, help='Matrix file with pairwise statistics between markers.')
	parser.add_option( '', '--out1', dest='out1', default='UPGMA_ordered_scaff.txt', help='Output file name of ordered scaffold. [default: %default]')
	parser.add_option( '', '--out2', dest='out2', default='UPGMA_ordered_mark.txt', help='Output file name of ordered markers. [default: %default]')
	parser.add_option( '', '--type', dest='type', default='IDENT', help='Data type : IDENT (marker linkage) or DIF (recombination rate). [default: %default]')
	(options, args) = parser.parse_args()
	
	if options.scaff == None:
		sys.exit('--scaff argument is missing')
	if options.mat == None:
		sys.exit('--mat argument is missing')
	
	#Recording mean distance/ressemblance between scaffold
	temp1 = tempfile.NamedTemporaryFile()
	mat_mark2mat_scaff(options.scaff, options.mat, temp1)
	temp1.flush()
	
	# Registering markers order
	#on enregistre dans une liste les markers dans l'ordre dans la matrice
	file = open(options.mat)
	index_mat_mark = file.readline().split()[1:]
	file.close()
	
	#Registering scaffold order and markers order
	#on enregistre dans un dico les scaffolds et l'ordre des positions (dans la matrice) des marqueurs dans le scaffold. Dans un second dictionnaire on enregistre l'ordre des scaffolds
	file = open(options.scaff)
	dico_mark = {}
	dico_orient = {}
	liste_ordre = []
	dico_index = {}
	nb_mark = 0
	for line in file:
		data = line.split()
		if data:
			nb_mark += 1
			if data[1] in dico_mark:
				dico_mark[data[1]].append(index_mat_mark.index(data[0]))
				dico_index[index_mat_mark.index(data[0])] = data[0]
			else:
				dico_mark[data[1]] = []
				dico_mark[data[1]].append(index_mat_mark.index(data[0]))
				dico_index[index_mat_mark.index(data[0])] = data[0]
				liste_ordre.append(data[1])
				dico_orient[data[1]] = 'FWD'
	file.close()
	nb_mark = float(nb_mark)
	
	#Recording markers pairwise information
	matrix_mark = {}
	info_scaff_mark = {}
	file = open(options.mat)
	file.readline()
	for line in file:
		data = line.split()
		if data:
			if data[0] in matrix_mark:
				sys.exit("Warning there is redundancy in markers names, verify your data")
			else:
				matrix_mark[data[0]] = map(float, data[1:])
	file.close()
	
	
	#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	#!!WARNING dico_orient is changing everytime in cherche_best_pos()!!
	#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	
	
	#Registering marker order in scaffold in a dictionnary
	file = open(temp1.name)
	index_mat = file.readline().split()[2:]
	file.close()
	
	#Recording scaffold pairwise information
	matrix = {}
	info_scaff = {}
	file = open(temp1.name)
	file.readline()
	for line in file:
		data = line.split()
		if data:
			if data[0] in matrix:
				sys.exit("Warning there is redundancy in markers names, verify your data")
			else:
				matrix[data[0]] = map(float, data[2:])
				info_scaff[data[0]] = int(data[1])
	file.close()
	
	#Doing the UPGMA
	GROUP = pseudoUP(matrix, index_mat, info_scaff, options.type, matrix_mark, liste_ordre, dico_mark, dico_orient, dico_index, nb_mark)
	
	#Outputing ordonned markers
	file = open(options.scaff)
	dico_mark = {}
	for line in file:
		data = line.split()
		if data:
			if data[1] in dico_mark:
				dico_mark[data[1]].append(data[0])
			else:
				dico_mark[data[1]] = []
				dico_mark[data[1]].append(data[0])
	outfile1 = open(options.out1,'w')
	outfile = open(options.out2,'w')
	for n in GROUP[1]:
		outfile1.write(n+'\t'+dico_orient[n]+'\n')
		liste = list(dico_mark[n])
		if dico_orient[n] == 'REV':
			liste.reverse()
		for j in liste:
			outfile.write(j+'\t'+n+'\n')
	outfile.close()
	outfile1.close()
	
if __name__ == "__main__": __main__()