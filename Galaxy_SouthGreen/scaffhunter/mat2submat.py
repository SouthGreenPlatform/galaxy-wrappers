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


def __main__():
	#Parse Command Line
	parser = optparse.OptionParser(usage="python %prog [options]\n\nProgram designed by Guillaume MARTIN : guillaume.martin@cirad.fr\n\n"
	"This program takes a square matrix (--mat) file and create de sub-matrix containing all IDs provided in another file (--mark).\n")
	# Wrapper options.
	parser.add_option( '', '--mark', dest='mark', default=None, help='A tabulated file containing in col 1 the id that will be contained in the sub-matrix')
	parser.add_option( '', '--mat', dest='mat', default=None, help='Square matrix file')
	parser.add_option( '', '--out', dest='out', default='sub_matrix.txt', help='Output file name. [default: %default]')
	(options, args) = parser.parse_args()
	
	if options.mark == None:
		sys.exit('--mark argument is missing')
	if options.mat == None:
		sys.exit('--mat argument is missing')
	
	#On enregistre les id de la matrice
	os.system('echo "Registering id matrix"')
	file = open(options.mat)
	index_mat = file.readline().split()[1:]
	file.close()
	
	#on enregistre dans un dico les position des markers a garder
	os.system('echo "Registering id to be kept"')
	file = open(options.mark)
	kept = []
	pos_kept = []
	for line in file:
		data = line.split()
		if data:
			if data[0] in kept:
				sys.exit("The programme exited without finishing : There is redundancy in markers names in the sub-id file")
			elif not(data[0] in index_mat):
				mot = "The programme exited without finishing : The marker "+data[0]+"is not in the matrix"
				sys.exit(mot)
			kept.append(data[0])
			pos_kept.append(index_mat.index(data[0]))
	file.close()
		
	#on charge la matrice dans un dico
	os.system('echo "Recording pairwise information"')
	matrix = {}
	file = open(options.mat)
	file.readline()
	for line in file:
		data = line.split()
		if data:
			if data[0] in matrix:
				sys.exit("The programme exited without finishing : There is redundancy in markers names in the matrix")
			else:
				matrix[data[0]] = data[1:]
	file.close()
	
	#Creating the table:
	os.system('echo "Creating the sub-table"')
	dico_table = {}
	for n in kept:
		dico_table[n] = []
		for j in pos_kept:
			dico_table[n].append(matrix[n][j])
		
	#Printing results
	os.system('echo "Printing results"')
	outfile = open(options.out,'w')
	outfile.write('ID\t'+'\t'.join(kept)+'\n')
	for n in kept:
		outfile.write(n+'\t'+'\t'.join(dico_table[n])+'\n')
	outfile.close()
	
if __name__ == "__main__": __main__()