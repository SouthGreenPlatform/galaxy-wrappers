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


def find_info(LINE):
	dic = {}
	for n in LINE[11:]:
		liste = n.split(':')
		if len(liste) > 2:
			dic[liste[0]] = liste[2]
	return dic

def CIGAR(mot):
	if 'N' in mot:
		sys.exit('N in mot')
	elif 'P' in mot:
		sys.exit('P in mot')
	elif '=' in mot:
		sys.exit('= in mot')
	elif 'X' in mot:
		sys.exit('X in mot')
	else:
		mot_split = mot.replace('M',' M ').replace('D',' D ').replace('I',' I ').replace('S',' S ').replace('H',' H ').split()
	total = 0
	while mot_split:
		if mot_split[1] == 'M' or mot_split[1] == 'D':
			total = total + int(mot_split[0])
		del mot_split[0]
		del mot_split[0]
	return total

def sam2tab(LINE):
	TOTAL = CIGAR(LINE[5])
	if LINE[1] == '0':
		return LINE[0]+'\t'+LINE[2]+'\t-\t1\t-\t'+LINE[3]+'\t'+str(int(LINE[3])+TOTAL)+'\t1'
	elif LINE[1] == '16':
		return LINE[0]+'\t'+LINE[2]+'\t-\t1\t-\t'+LINE[3]+'\t'+str(int(LINE[3])+TOTAL)+'\t-1'
	else:
		sys.exit('Probleme in sam flag '+LINE[1])

def __main__():
	#Parse Command Line
	parser = optparse.OptionParser(usage="python %prog [options]\n\nProgram designed by Guillaume MARTIN : guillaume.martin@cirad.fr")
	# Wrapper options. 
	parser.add_option( '', '--sam', dest='sam', help='The sam file')
	parser.add_option( '', '--dif', dest='dif', help='The minimal difference between the best and second hit accepted to consider the hit as single')
	parser.add_option( '', '--type', dest='type', help='The output type, sam or tab')
	(options, args) = parser.parse_args()
	
	min_dif = int(options.dif)
	
	file = open(options.sam)
	if options.type != 'sam' and options.type != 'tab':
		sys.exit('Output option is not recognised')
	for line in file:
		data = line.split()
		if data:
			if data[0][0] != '@':
				if data[1] != '256' and data[1] != '272':
					dico = find_info(data)
					if 'XS' in dico:
						if abs(int(dico['AS'])-int(dico['XS'])) >= min_dif:
							if options.type == 'sam':
								print('\t'.join(data))
							else:
								if data[1] != '4':
									print(sam2tab(data))
					else:
						if options.type == 'sam':
							print('\t'.join(data))
						else:
							if data[1] != '4':
								print(sam2tab(data))
			else:
				if options.type == 'sam':
					print('\t'.join(data))
	
if __name__ == "__main__": __main__()



