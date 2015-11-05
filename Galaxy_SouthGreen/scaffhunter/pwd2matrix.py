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
	"This program takes a file containing pairwise statistics and returns a matrix containing the same information.")
	# Wrapper options. 
	parser.add_option( '', '--pwd', dest='pwd', default=None, help='Pairwise file. Column 1: id1, column 2: id2, column 3: statistics')
	parser.add_option( '', '--out', dest='out', default='matrix.tab', help='Output file name. [default: %default]')
	
	(options, args) = parser.parse_args()
	
	if options.pwd == None:
		sys.exit('Please provide a pairwise file in --pwd argument')
	
	#recording markers ids
	os.system('echo "recording markers ids"')
	liste_id = []
	file = open(options.pwd)
	for line in file:
		data=line.split()
		if data:
			if not(data[0] in liste_id):
				liste_id.append(data[0])
			if not(data[1] in liste_id):
				liste_id.append(data[1])
	file.close()
	
	#creating the table structure
	os.system('echo "creating the table structure"')
	liste_vide = []
	for n in liste_id:
		liste_vide.append('999999999')
	dic_LOD = {}
	for n in liste_id:
		dic_LOD[n] = list(liste_vide)
	
	#filling the table
	os.system('echo "filling the table"')
	file = open(options.pwd)
	for line in file:
		data = line.split()
		if data:
			dic_LOD[data[0]][liste_id.index(data[1])] = float(data[2])
			dic_LOD[data[1]][liste_id.index(data[0])] = float(data[2])
	file.close()
				
	#writing results
	os.system('echo "writing results"')
	outfile1 = open(options.out,'w')
	outfile1.write('ID\t'+'\t'.join(liste_id)+'\n')
	for n in liste_id:
		outfile1.write(n+'\t'+'\t'.join(map(str,dic_LOD[n]))+'\n')
	outfile1.close()
	
if __name__ == "__main__": __main__()