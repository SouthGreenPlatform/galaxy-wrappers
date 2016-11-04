
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

def factorielle(x):
	if x < 2:
		return 1
	else:
		return x * factorielle(x-1)

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

def rear_fait(ORD, ORIENT):
	mot = ''.join(ORD)
	for n in ORD:
		mot = mot+ORIENT[n]
	return mot

def reorderient(ORDRE, ORIENT):
	NEW_ORDRE = list(ORDRE)
	NEW_ORIENT = {}
	for n in ORIENT:
		NEW_ORIENT[n] = ORIENT[n]
	if random.randint(0,1) == 0:
		#On travail sur 1 scaffold
		pos = random.randint(0, (len(ORDRE)-1))
		scaffold = NEW_ORDRE[pos]
		if random.randint(0,1) == 0:
			#On fait un deplacement de scaffold
			pos_new = pos
			while pos_new == pos:
				pos_new = random.randint(0, (len(ORDRE)-1))
			del NEW_ORDRE[pos]
			if scaffold in NEW_ORDRE:
				sys.exit('bug')
			NEW_ORDRE.insert(pos_new, scaffold)
			if random.randint(0,1) == 0:
				#On fait en plus une reorientation
				if NEW_ORIENT[scaffold] == 'FWD':
					NEW_ORIENT[scaffold] = 'REV'
				elif NEW_ORIENT[scaffold] == 'REV':
					NEW_ORIENT[scaffold] = 'FWD'
				else:
					sys.exit('bug')
		else:
			#On fait juste une reorientation
			if NEW_ORIENT[scaffold] == 'FWD':
				NEW_ORIENT[scaffold] = 'REV'
			elif NEW_ORIENT[scaffold] == 'REV':
				NEW_ORIENT[scaffold] = 'FWD'
			else:
				sys.exit('bug')
	else:
		pos_debut = random.randint(0, (len(ORDRE)-2))
		pos_fin = pos_debut
		while (pos_debut + 1) >= pos_fin and len(NEW_ORDRE[pos_debut:pos_fin]) == len(NEW_ORDRE):
			pos_fin = random.randint(1, len(ORDRE))
		scaffold = NEW_ORDRE[pos_debut:pos_fin]
		del NEW_ORDRE[pos_debut:pos_fin]
		if random.randint(0,1) == 0:
			#On fait un deplacement de scaffold
			pos_new = pos_debut
			while pos_new == pos_debut:
				pos_new = random.randint(0, len(NEW_ORDRE))
			if random.randint(0,1) == 0:
				#On fait en plus une reorientation
				for n in scaffold:
					NEW_ORDRE.insert(pos_new, n)
					if NEW_ORIENT[n] == 'FWD':
						NEW_ORIENT[n] = 'REV'
					elif NEW_ORIENT[n] == 'REV':
						NEW_ORIENT[n] = 'FWD'
					else:
						sys.exit('bug')
			else:
				#On fait juste un deplacement
				for n in scaffold:
					NEW_ORDRE.insert(pos_new, n)
					pos_new += 1
		else:
			#On fait juste une reorientation
			for n in scaffold:
				NEW_ORDRE.insert(pos_debut, n)
				if NEW_ORIENT[n] == 'FWD':
					NEW_ORIENT[n] = 'REV'
				elif NEW_ORIENT[n] == 'REV':
					NEW_ORIENT[n] = 'FWD'
				else:
					sys.exit('bug')
		if len(NEW_ORDRE) != len(ORDRE):
			sys.exit('bug')
	return [NEW_ORDRE, NEW_ORIENT]

def __main__():
	#Parse Command Line
	parser = optparse.OptionParser(usage="python %prog [options]\n\nProgram designed by Guillaume MARTIN : guillaume.martin@cirad.fr\n\n"
	"This program take in input a matrix containing maker linkage or divergence and the two output files of UPGMA.py and try to optimize the "
	"scaffold order and orientation proposed by UPGMA.py.\n"
	"The optimization is performed by calculating a score for the scaffold order, trying rearrangements (scaffold or scaffold group permutations "
	"and/or inversion) followed by score re-calculation. If the new score is better than the previous, the new order is conserved; else, the previous "
	"order is conserved. The program stops when more than a defined number of consecutive rearrangement (passed in --iter argument) do not improve the ordering.\n"
	"This program output a file containing ordered scaffolds in a table format and a table file containing ordered markers relative to the scaffold order.\n"
	)
	# Wrapper options.
	parser.add_option( '', '--mat', dest='mat', default=None, help='Matrix file with pairwise statistics between markers.')
	parser.add_option( '', '--scaff', dest='scaff', default=None, help='The output file passed in --out1 argument of UPGMA.py. This file contains ordered and orientated scaffolds ate the end of UPGMA.py program.')
	parser.add_option( '', '--mark', dest='mark', default=None, help='A table file containing in column 1: markers name, column2: scaffold name. This is the file generated by UPGMA.py in --out2 argument')
	parser.add_option( '', '--iter', dest='iter', default='auto', help='Number of consecutive rearrangement tried without improvement to stop the improvement of the ordering and orientation (integer or auto). [default: %default]')
	parser.add_option( '', '--type', dest='type', default='IDENT', help='Data type : IDENT (marker linkage) or DIF (recombination rate). [default: %default]')
	parser.add_option( '', '--out1', dest='out1', default='UPGMA_ordered_opt_scaf.txt', help='Output file name of ordered scaffold. [default: %default]')
	parser.add_option( '', '--out2', dest='out2', default='UPGMA_ordered_opt_mark.txt', help='Output file name of ordered markers. [default: %default]')
	(options, args) = parser.parse_args()
	
	if options.scaff == None:
		sys.exit('--scaff argument is missing')
	if options.mark == None:
		sys.exit('--mark argument is missing')
	if options.mat == None:
		sys.exit('--mat argument is missing')
	
	
	#on enregistre dans une liste les markers dans l'ordre dans la matrice
	os.system('echo "Registering scaffold order"')
	file = open(options.mat)
	index_mat = file.readline().split()[1:]
	file.close()
	
	# On enregistre l'orientation des scaffolds
	os.system('echo "Registering scaffold orientation"')
	dico_orient = {}
	file = open(options.scaff)
	for line in file:
		data = line.split()
		if data:
			dico_orient[data[0]] = data[1]
	file.close()
	
	#on enregistre dans un dico les scaffolds et l'ordre des positions (dans la matrice) des marqueurs dans le scaffold. Dans une liste on enregistre l'ordre des scaffolds
	os.system('echo "Registering scaffold order and markers order"')
	file = open(options.mark)
	dico_mark = {}#pour chaque cle = scaffold, on enregistre l'ordre des markers dans ce scaffold
	liste_ordre = []#Ordre des scaffolds
	dico_index = {}#position du markers dans la matrice
	# dico_orient = {}#######
	nb_mark = 0
	for line in file:
		data = line.split()
		if data:
			nb_mark += 1
			if data[1] in dico_mark:
				dico_index[index_mat.index(data[0])] = data[0]
				if dico_orient[data[1]] == 'FWD':
					dico_mark[data[1]].append(index_mat.index(data[0]))
				elif dico_orient[data[1]] == 'REV':
					dico_mark[data[1]].insert(0, index_mat.index(data[0]))
				else:
					sys.exit('bugggggg')
			else:
				dico_mark[data[1]] = []
				dico_mark[data[1]].append(index_mat.index(data[0]))
				dico_index[index_mat.index(data[0])] = data[0]
				liste_ordre.append(data[1])
				# dico_orient[data[1]] = 'FWD'################
	file.close()
	nb_mark = float(nb_mark)
	
	#On verifie le nombre d'iteration
	total_possible = ((float(len(liste_ordre)-1)**2)*2)+(len(liste_ordre))
	if options.iter == 'auto':
		iteration = total_possible
		os.system('echo "The number of iteration is automatically estimated to the number of possible iterations in case of one scaffold repositioned : '+str(total_possible)+'"')
	elif options.iter == 'all':
		os.system('echo "All possible rearrangments are tested"')
	else:
		iteration = int(options.iter)
		if iteration > total_possible:
			os.system('echo "The number of iteration is over the number of possible iterations in case of one scaffold repositioned. It has been ajusted to : '+str(total_possible)+' if not the script can be in a infinite loop"')
			iteration = total_possible
		else:
			os.system('echo "The number of iteration is under the number of possible iterations in case of one scaffold repositioned ('+str(total_possible)+'). You choose : '+str(iteration)+'"')
	
	#on charge la matrice dans un dico
	os.system('echo "Recording pairwise information"')
	matrix = {}
	info_scaff = {}
	file = open(options.mat)
	file.readline()
	for line in file:
		data = line.split()
		if data:
			if data[0] in matrix:
				sys.exit("Warning there is redundancy in markers names, verify your data")
			else:
				matrix[data[0]] = map(float, data[1:])
	file.close()
	
	if options.type == 'IDENT':
		#On essaye toutes les orientations possibles de cet ordre
		os.system('echo "Looking for best orientation in this order"')
		t0 = time.clock()
		score = calcul_score_id(matrix, liste_ordre, dico_mark, dico_orient, dico_index, nb_mark)
		os.system('echo "Reference score : '+str(score)+'; time spent in calculation : '+str(time.clock()-t0)+'"')
		for n in liste_ordre:
			t0 = time.clock()
			dico_reor = {}
			for j in dico_orient:
				if j == n:
					if dico_orient[j] == 'FWD':
						dico_reor[j] = 'REV'
					elif dico_orient[j] == 'REV':
						dico_reor[j] = 'FWD'
					else:
						sys.exit('bug')
					# print n, dico_reor[j], dico_orient[n]
				else:
					dico_reor[j] = str(dico_orient[j])
			score_reor = calcul_score_id(matrix, liste_ordre, dico_mark, dico_reor, dico_index, nb_mark)
			os.system('echo "Reorient score : '+str(score_reor)+'; time spent in calculation : '+str(time.clock()-t0)+'"')
			if score < score_reor:
				dico_orient = {}
				for k in dico_reor:
					dico_orient[k] = str(dico_reor[k])
				score = score_reor
		#On fait des rearrangements
		os.system('echo "Performing rearrangements"')
		iter = 0
		fait = set()
		fait.add(rear_fait(liste_ordre, dico_orient))
		while iter < iteration:
			t0 = time.clock()
			iter += 1
			REOR = reorderient(liste_ordre,dico_orient)
			j = 0
			while rear_fait(REOR[0], REOR[1]) in fait:
				j += 1
				REOR = reorderient(liste_ordre,dico_orient)
			fait.add(rear_fait(REOR[0], REOR[1]))
			score_reor = calcul_score_id(matrix, REOR[0], dico_mark, REOR[1], dico_index, nb_mark)
			os.system('echo "iteration : '+str(iter)+'; Reoriented score : '+str(score_reor)+', best score : '+str(score)+'; time spent in calculation : '+str(time.clock()-t0)+'"')
			if score < score_reor:
				dico_orient = {}
				for n in REOR[1]:
					dico_orient[n] = REOR[1][n]
				score = score_reor
				liste_ordre = list(REOR[0])
				iter = 0
				fait = set()
				fait.add(rear_fait(liste_ordre, dico_orient))
		#On fait une derniere etape pour voir ce qui est oriente et ce qui ne l'est pas
		os.system('echo "Looking for orientated scaffolds and those who are not"')
		os.system('echo "Reference score : '+str(score)+'; time spent in calculation : '+str(time.clock()-t0)+'"')
		dico_ordre_final = {}
		for n in liste_ordre:
			t0 = time.clock()
			dico_reor = {}
			nom_ord = n
			for j in dico_orient:
				if j == n:
					if dico_orient[j] == 'FWD':
						dico_reor[j] = 'REV'
					else:
						dico_reor[j] = 'FWD'
				else:
					dico_reor[j] = dico_orient[j]
			score_reor = calcul_score_id(matrix, liste_ordre, dico_mark, dico_reor, dico_index, nb_mark)
			os.system('echo "Reorient score : '+str(score_reor)+'; time spent in calculation : '+str(time.clock()-t0)+'"')
			if score < score_reor:
				dico_orient = {}
				for k in dico_reor:
					dico_orient[k] = dico_reor[k]
				score = score_reor
				dico_ordre_final[n] = 'Ord'
			elif score == score_reor:
				dico_ordre_final[n] = 'NoOrd'
			else:
				dico_ordre_final[n] = 'Ord'
	

	elif options.type == 'DIF':
		#On essaye toutes les orientations possibles de cet ordre
		os.system('echo "Looking for best orientation in this order"')
		t0 = time.clock()
		score = calcul_score_dif(matrix, liste_ordre, dico_mark, dico_orient, dico_index, nb_mark)
		os.system('echo "Reference score : '+str(score)+'; time spent in calculation : '+str(time.clock()-t0)+'"')
		for n in liste_ordre:
			t0 = time.clock()
			dico_reor = {}
			for j in dico_orient:
				if j == n:
					if dico_orient[j] == 'FWD':
						dico_reor[j] = 'REV'
					elif dico_orient[j] == 'REV':
						dico_reor[j] = 'FWD'
					else:
						sys.exit('bug')
				else:
					dico_reor[j] = dico_orient[j]
			score_reor = calcul_score_dif(matrix, liste_ordre, dico_mark, dico_reor, dico_index, nb_mark)
			os.system('echo "Reorient score : '+str(score_reor)+'; time spent in calculation : '+str(time.clock()-t0)+'"')
			if score < score_reor:
				dico_orient = {}
				for k in dico_reor:
					dico_orient[k] = dico_reor[k]
				score = score_reor
		#On fait des rearrangements
		os.system('echo "Performing rearrangements"')
		iter = 0
		fait = set()
		fait.add(rear_fait(liste_ordre, dico_orient))
		while iter < iteration:
			t0 = time.clock()
			iter += 1
			REOR = reorderient(liste_ordre,dico_orient)
			while rear_fait(REOR[0], REOR[1]) in fait:
				REOR = reorderient(liste_ordre,dico_orient)
			fait.add(rear_fait(REOR[0], REOR[1]))
			score_reor = calcul_score_dif(matrix, REOR[0], dico_mark, REOR[1], dico_index, nb_mark)
			os.system('echo "iteration : '+str(iter)+'; Reoriented score : '+str(score_reor)+', best score : '+str(score)+'; time spent in calculation : '+str(time.clock()-t0)+'"')
			if score < score_reor:
				dico_orient = {}
				for n in REOR[1]:
					dico_orient[n] = REOR[1][n]
				score = score_reor
				liste_ordre = list(REOR[0])
				iter = 0
				fait = set()
				fait.add(rear_fait(liste_ordre, dico_orient))
		#On fait une derniere etape pour voir ce qui est oriente et ce qui ne l'est pas
		os.system('echo "Looking for orientated scaffolds and those who are not"')
		os.system('echo "Reference score : '+str(score)+'; time spent in calculation : '+str(time.clock()-t0)+'"')
		dico_ordre_final = {}
		for n in liste_ordre:
			t0 = time.clock()
			dico_reor = {}
			nom_ord = n
			for j in dico_orient:
				if j == n:
					if dico_orient[j] == 'FWD':
						dico_reor[j] = 'REV'
					else:
						dico_reor[j] = 'FWD'
				else:
					dico_reor[j] = dico_orient[j]
			score_reor = calcul_score_dif(matrix, liste_ordre, dico_mark, dico_reor, dico_index, nb_mark)
			os.system('echo "Reorient score : '+str(score_reor)+'; time spent in calculation : '+str(time.clock()-t0)+'"')
			if score < score_reor:
				dico_orient = {}
				for k in dico_reor:
					dico_orient[k] = dico_reor[k]
				score = score_reor
				dico_ordre_final[n] = 'Ord'
			elif score == score_reor:
				dico_ordre_final[n] = 'NoOrd'
				if dico_orient[n] == 'REV':				#
					dico_orient = {}					#
					for k in dico_reor:					#
						dico_orient[k] = dico_reor[k]	#
			else:
				dico_ordre_final[n] = 'Ord'
	else:
		sys.exit('Wrong argument passed in --type')
	
	#On genere le fichier pour matrix to orthodother
	os.system('echo "Outputing ordonned markers"')
	outfile1 = open(options.out1,'w')
	outfile = open(options.out2,'w')
	for n in liste_ordre:
		outfile1.write(n+'\t'+dico_orient[n]+'\t'+dico_ordre_final[n]+'\n')
		liste = list(dico_mark[n])
		if dico_orient[n] == 'REV':
			liste.reverse()
		for j in liste:
			outfile.write(dico_index[j]+'\t'+n+'\n')
	outfile.close()
	outfile1.close()
	
if __name__ == "__main__": __main__()