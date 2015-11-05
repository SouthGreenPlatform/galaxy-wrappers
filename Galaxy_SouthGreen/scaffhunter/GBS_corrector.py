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

import optparse, os, shutil, subprocess, sys, tempfile, fileinput, ConfigParser, operator, time


def cherche_best(DICO):
	BEST = 'Un'
	num = 0
	for n in DICO:
		if DICO[n] > num:
			BEST = n
			num = DICO[n]
		elif DICO[n] == num:
			BEST = 'Un'
	return BEST

def correct(HEADER, LISTE_MARK, DICO_IND_ERROR, DICO_MARK_ERROR, DICO, FENETRE):
	#We correct by individuals
	for ind in HEADER:
		for mark in LISTE_MARK:
			if LISTE_MARK.index(mark) < FENETRE:
				dico_var = {}
				for n in LISTE_MARK[0:(2*FENETRE)+1]:
					if DICO[n][HEADER.index(ind)] in dico_var:
						dico_var[DICO[n][HEADER.index(ind)]] = dico_var[DICO[n][HEADER.index(ind)]] + 1
					else:
						dico_var[DICO[n][HEADER.index(ind)]] = 1
				best = cherche_best(dico_var)
				if best != DICO[mark][HEADER.index(ind)] and best != 'Un' and best != '--':
					DICO_MARK_ERROR[mark] = DICO_MARK_ERROR[mark] + 1
					DICO_IND_ERROR[ind] = DICO_IND_ERROR[ind] + 1
					DICO[mark][HEADER.index(ind)] = best
			elif LISTE_MARK.index(mark) >= (len(LISTE_MARK) - FENETRE):
				dico_var = {}
				for n in LISTE_MARK[((len(LISTE_MARK)-2*FENETRE)-1):len(LISTE_MARK)]:
					if DICO[n][HEADER.index(ind)] in dico_var:
						dico_var[DICO[n][HEADER.index(ind)]] = dico_var[DICO[n][HEADER.index(ind)]] + 1
					else:
						dico_var[DICO[n][HEADER.index(ind)]] = 1
				best = cherche_best(dico_var)
				if best != DICO[mark][HEADER.index(ind)] and best != 'Un' and best != '--':
					DICO_MARK_ERROR[mark] = DICO_MARK_ERROR[mark] + 1
					DICO_IND_ERROR[ind] = DICO_IND_ERROR[ind] + 1
					DICO[mark][HEADER.index(ind)] = best
			else:
				dico_var = {}
				for n in LISTE_MARK[LISTE_MARK.index(mark)-FENETRE:LISTE_MARK.index(mark)+FENETRE+1]:
					if DICO[n][HEADER.index(ind)] in dico_var:
						dico_var[DICO[n][HEADER.index(ind)]] = dico_var[DICO[n][HEADER.index(ind)]] + 1
					else:
						dico_var[DICO[n][HEADER.index(ind)]] = 1
				best = cherche_best(dico_var)
				if best != DICO[mark][HEADER.index(ind)] and best != 'Un' and best != '--':
					DICO_MARK_ERROR[mark] = DICO_MARK_ERROR[mark] + 1
					DICO_IND_ERROR[ind] = DICO_IND_ERROR[ind] + 1
					DICO[mark][HEADER.index(ind)] = best
	return [DICO_IND_ERROR, DICO_MARK_ERROR, DICO]
	

def __main__():
	#Parse Command Line
	parser = optparse.OptionParser(usage="python %prog [options]\n\nProgram designed by Guillaume MARTIN : guillaume.martin@cirad.fr\n\n"
	"This program identifies markers genotyping errors recorded in a table file based on their order (obtained from genetic map or reference "
	"sequence) provided in a table file. This program is based on the principle that no more than one recombination event can occur in a window "
	"of x around the observed marker (x is passed in --fen argument). If more than one recombination event is observed, a genotyping error is "
	"identified. This program output three files:\n"
	"\t1) A file containing corrected markers (--nosu argument)\n"
	"\t2) A file containing corrected markers where redundancy is removed (--nr argument)\n"
	"\t3) A file where markers with more than X errors are identified are removed. (--Nosu argument)\n"
	"This program works on phased markers of same type. It is based on the principle that no more than one recombination can be observed in a window around each markers.")
	
	# Wrapper options. 
	parser.add_option( '', '--table', dest='table', default=None, help='The table file containing phased data. col 1 : markers name, col 2, 3, 4 necessary but not used by the program, col 5 to end : individual genotypes. First line contain header. Redundant names (markers and individuals) are not allowed')
	parser.add_option( '', '--fen', dest='fen', default="10=3", help='Minimal marker number in a linkage group/scaffold/reference sequence to apply correction followed by half window size. Both values should be separated by "=". If different classes, couple can be repeated and separated by "=". [default: %default]')
	parser.add_option( '', '--order', dest='order', default=None, help='A table file with markers names in column 1 and scaffold they match with in column2. Markers should be ordered first on linkage group/scaffold/reference sequence and second on markers position on linkage group/scaffold/reference sequence')
	parser.add_option( '', '--suspect', dest='suspect', default=10, help='The minimal number of correction in a marker to consider this marker as supect. [default: %default]')
	parser.add_option( '', '--nr', dest='nr', default='corrected_non_redunddant_mark.tab', help='The name of the corrected non redundant output file. [default: %default]')
	parser.add_option( '', '--cor', dest='cor', default='corrected_mark.tab', help='The name of the corrected output file. [default: %default]')
	parser.add_option( '', '--Nosu', dest='Nosu', default='filtered_mark.tab', help='The name of the output file containing the non suspect markers. [default: %default]')
	(options, args) = parser.parse_args()
	
	
	if options.table == None:
		sys.exit('--table argument is missing')
	if options.order == None:
		sys.exit('--order argument is missing')
	
	susp = int(options.suspect)
	
	fenetre = []
	toto = options.fen.split("=")
	while toto:
		fenetre.append([int(toto[0]),int(toto[1])])
		del toto[0]
		del toto[0]
	
	fenetre.sort(key=operator.itemgetter(0), reverse=True)
	# print fenetre
	
	#recording data in dictionnary
	dico = {}
	dico_3_col = {}
	dico_mark_error = {}
	dico_ind_error = {}
	header = []
	file = open(options.table)
	line = file.readline()
	header = line.split()
	taille_line = len(header)
	header = line.split()[4:]
	head_3_col = line.split()[1:4]
	liste_mark = []
	print ('Your file contain '+str(taille_line-4)+' individuals')
	for line in file:
		data = line.split()
		if data:
			if len(data) != taille_line:
				sys.exit('Please check your --table file, lines have different length')
			if data[0] in liste_mark:
				sys.exit('Please check your --table file, there is redundancy in marker names')
			dico[data[0]] = data[4:]
			dico_3_col[data[0]] = data[1:4]
			liste_mark.append(data[0])
			dico_mark_error[data[0]] = 0
	for ind in header:
		dico_ind_error[ind] = 0
	
	#On travail maintenant par scaffold
	liste_mark = []
	sub_liste_mark = []
	rev_dic_mark = {}
	file = open(options.order)
	scaffold = ''
	for line in file:
		data = line.split()
		if data:
			if data[0] in liste_mark:
				sys.exit('Please check your --order file, there is redundancy in marker names')
			if scaffold == '':
				liste_mark.append(data[0])
				sub_liste_mark.append(data[0])
				scaffold = data[1]
				rev_dic_mark[data[0]] = data[1]
			elif scaffold == data[1]:
				liste_mark.append(data[0])
				sub_liste_mark.append(data[0])
				rev_dic_mark[data[0]] = data[1]
			else:
				for k in fenetre:
					if k[0] <= len(sub_liste_mark):
						print scaffold, len(sub_liste_mark), k
						resultat = correct(header, sub_liste_mark, dico_ind_error, dico_mark_error, dico, k[1])
						dico_ind_error = resultat[0]
						dico_mark_error = resultat[1]
						dico = resultat[2]
						break
				sub_liste_mark = []
				sub_liste_mark.append(data[0])
				liste_mark.append(data[0])
				scaffold = data[1]
				rev_dic_mark[data[0]] = data[1]
	if sub_liste_mark:
		for k in fenetre:
			if k[0] <= len(sub_liste_mark):
				print scaffold, len(sub_liste_mark), k
				resultat = correct(header, sub_liste_mark, dico_ind_error, dico_mark_error, dico, k[1])
				dico_ind_error = resultat[0]
				dico_mark_error = resultat[1]
				dico = resultat[2]
				break
				
	#now we print output
	outfile_nr = open(options.nr,'w')
	outfile = open(options.cor,'w')
	outfile.write('marker')
	outfile_nr.write('marker')
	for n in head_3_col:
		outfile.write('\t'+n)
		outfile_nr.write('\t'+n)
	outfile.write('\tmark_error')
	outfile_nr.write('\tmark_error')
	for n in header:
		outfile.write('\t'+n)
		outfile_nr.write('\t'+n)
	outfile.write('\n')
	outfile.write('ind_error\t--\t-\t-\t-')
	outfile_nr.write('\n')
	outfile_nr.write('ind_error\t--\t-\t-\t-')
	for n in header:
		outfile.write('\t'+str(dico_ind_error[n]))
		outfile_nr.write('\t'+str(dico_ind_error[n]))
	outfile.write('\n')
	outfile_nr.write('\n')
	
	mark_prec = []
	scaf_prec = ''
	for mark in liste_mark:
		outfile.write(mark)
		if dico[mark] == mark_prec and rev_dic_mark[mark] == scaf_prec:
			mark_prec = dico[mark]
			scaf_prec = rev_dic_mark[mark]
		else:
			outfile_nr.write(mark)
			mark_prec = dico[mark]
			scaf_prec = rev_dic_mark[mark]
			for n in dico_3_col[mark]:
				outfile_nr.write('\t'+n)
			outfile_nr.write('\t'+str(dico_mark_error[mark]))
			for ind in header:
				outfile_nr.write('\t'+dico[mark][header.index(ind)])
			outfile_nr.write('\n')
		for n in dico_3_col[mark]:
			outfile.write('\t'+n)
		outfile.write('\t'+str(dico_mark_error[mark]))
		for ind in header:
			outfile.write('\t'+dico[mark][header.index(ind)])
		outfile.write('\n')
	outfile.close()
	outfile_nr.close()
	
	#Now we print the original file with suspect markers removed
	outfile = open(options.Nosu,'w')
	file = open(options.table)
	line = file.readline()
	outfile.write(line)
	for line in file:
		data = line.split()
		if data:
			if dico_mark_error[data[0]] < susp and data[0] in liste_mark:
				outfile.write(line)
	outfile.close()
	
if __name__ == "__main__": __main__()