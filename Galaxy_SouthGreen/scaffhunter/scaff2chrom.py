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


import optparse, os, shutil, subprocess, sys, tempfile, fileinput, ConfigParser, operator
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def rev_seq(seq):
	#function that reverse and complement a sequence
	my_dna = Seq(seq, generic_dna)
	return str(my_dna.reverse_complement())


def __main__():
	#Parse Command Line
	parser = optparse.OptionParser(usage="python %prog [options]\n\nProgram designed by Guillaume MARTIN : guillaume.martin@cirad.fr\n\n"
	"This program takes in input files a multifasta containing scaffolds and a table file containing scaffolds order calculated by reorderient.py and "
	"output a multifasta file containing reconstructed chromosomes and an agp file locating scaffolds into chromosomes.\n")
	
	# Wrapper options. 
	parser.add_option( '', '--table', dest='table', default=None, help='Table file containing in column 1: chromosome name, column 2: scaffold name, column 3: expected orientation (FWD or REV).')
	parser.add_option( '', '--seq', dest='seq', default=None, help='Multifasta sequence file containing scaffolds.')
	parser.add_option( '', '--unknown', dest='unknown', default='yes', help='Build an unknown chromosome with the remaining sequences: yes or no, [default: %default].')
	parser.add_option( '', '--out', dest='out', default='chromosomes.fasta', help='Fasta output file name.')
	parser.add_option( '', '--agp', dest='agp', default='chromosomes.agp', help='agp output file name.')
	(options, args) = parser.parse_args()
	
	if options.table == None:
		sys.exit('--table argument is missing')
	if options.seq == None:
		sys.exit('--seq argument is missing')
		
	#loading sequences
	record_dict = SeqIO.index(options.seq, "fasta")
	file = open(options.table)
	dico = set()
	chr = ''
	debut = 0
	sequence = ''
	scaff = ''
	outfile = open(options.out,'w')
	outfile2 = open(options.agp,'w')
	outfile2.write("##agp-version   2.0\n")
	status = ""
	for line in file:
		data = line.split()
		if data:
			if data[1] in dico:
				print 'Warning '+data[1]+' is found more than once'
			dico.add(data[1])
			if not(data[1] in record_dict):
				print data[1]+' is not in the fasta file'
				i -= 1
			elif chr == '':
				i = 0
				chr = data[0]
				if len(data) == 2:
					sequence = str(record_dict[data[1]].seq)
				else:#contain information on orientation
					if data[2] == 'FWD':
						sequence = str(record_dict[data[1]].seq)
					elif data[2] == 'REV':
						sequence = rev_seq(str(record_dict[data[1]].seq))
					else:
						sys.exit('The orientation information is not recognized')
				scaff = data[1]
				if len(data) >= 3:
					status = data[2]
				else:
					status = ""
			elif chr == data[0]:
				i += 1
				if status == 'FWD':
					outfile2.write(chr+'\t'+str(debut+1)+'\t'+str(len(sequence))+'\t'+str(i)+'\tW\t'+scaff+'\t1\t'+str(len(sequence)-debut)+'\t+\n')
				elif status == 'REV':
					outfile2.write(chr+'\t'+str(debut+1)+'\t'+str(len(sequence))+'\t'+str(i)+'\tW\t'+scaff+'\t1\t'+str(len(sequence)-debut)+'\t-\n')
				else:
					outfile2.write(chr+'\t'+str(debut+1)+'\t'+str(len(sequence))+'\t'+str(i)+'\tW\t'+scaff+'\t1\t'+str(len(sequence)-debut)+'\t+\n')
				debut = len(sequence)
				i += 1
				sequence = sequence + 'NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN'
				outfile2.write(chr+'\t'+str(debut+1)+'\t'+str(len(sequence))+'\t'+str(i)+'\tN\t100\tfragment\tno\n')
				debut = len(sequence)
				if len(data) == 2:
					sequence = sequence + str(record_dict[data[1]].seq)
				else:#contain information on orientation
					if data[2] == 'FWD':
						sequence = sequence + str(record_dict[data[1]].seq)
					elif data[2] == 'REV':
						sequence = sequence + rev_seq(str(record_dict[data[1]].seq))
					else:
						sys.exit('The orientation information is not recognized')
				scaff = data[1]
				if len(data) >= 3:
					status = data[2]
				else:
					status = ""
			else:
				i += 1
				if status == 'FWD':
					outfile2.write(chr+'\t'+str(debut+1)+'\t'+str(len(sequence))+'\t'+str(i)+'\tW\t'+scaff+'\t1\t'+str(len(sequence)-debut)+'\t+\n')
				elif status == 'REV':
					outfile2.write(chr+'\t'+str(debut+1)+'\t'+str(len(sequence))+'\t'+str(i)+'\tW\t'+scaff+'\t1\t'+str(len(sequence)-debut)+'\t-\n')
				else:
					outfile2.write(chr+'\t'+str(debut+1)+'\t'+str(len(sequence))+'\t'+str(i)+'\tW\t'+scaff+'\t1\t'+str(len(sequence)-debut)+'\t+\n')
				SeqIO.write(SeqRecord(Seq(sequence, generic_dna), id = chr, description=''),outfile, "fasta")
				chr = data[0]
				debut = 0
				i = 0
				sequence = ''
				if len(data) == 2:
					sequence = sequence + str(record_dict[data[1]].seq)
				else:#contain information on orientation
					if data[2] == 'FWD':
						sequence = sequence + str(record_dict[data[1]].seq)
					elif data[2] == 'REV':
						sequence = sequence + rev_seq(str(record_dict[data[1]].seq))
					else:
						sys.exit('The orientation information is not recognized')
				scaff = data[1]
				if len(data) >= 3:
					status = data[2]
				else:
					status = ""
	i += 1
	if status == 'FWD':
		outfile2.write(chr+'\t'+str(debut+1)+'\t'+str(len(sequence))+'\t'+str(i)+'\tW\t'+scaff+'\t1\t'+str(len(sequence)-debut)+'\t+\n')
	elif status == 'REV':
		outfile2.write(chr+'\t'+str(debut+1)+'\t'+str(len(sequence))+'\t'+str(i)+'\tW\t'+scaff+'\t1\t'+str(len(sequence)-debut)+'\t-\n')
	else:
		outfile2.write(chr+'\t'+str(debut+1)+'\t'+str(len(sequence))+'\t'+str(i)+'\tW\t'+scaff+'\t1\t'+str(len(sequence)-debut)+'\t+\n')
	SeqIO.write(SeqRecord(Seq(sequence, generic_dna), id = chr, description=''),outfile, "fasta")
	#all chromosomes have been constructed, now it is time to build unknown chromosome
	if options.unknown == 'yes':
		liste = []
		for n in record_dict:
			if not(n in dico):
				liste.append([n, len(str(record_dict[n].seq))])
		liste = sorted(liste, key=operator.itemgetter(1), reverse = True)
		sequence = ''
		debut = 0
		i = 0
		for n in liste:
			if sequence == '':
				sequence = str(record_dict[n[0]].seq)
			else:
				i += 1
				outfile2.write('chrUn_random\t'+str(debut+1)+'\t'+str(len(sequence))+'\t'+str(i)+'\tW\t'+scaff+'\t1\t'+str(len(sequence)-debut)+'\t+\n')
				debut = len(sequence)
				i += 1
				sequence = sequence + 'NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN'
				outfile2.write('chrUn_random\t'+str(debut+1)+'\t'+str(len(sequence))+'\t'+str(i)+'\tN\t100\tfragment\tno\n')
				debut = len(sequence)
				sequence = sequence + str(record_dict[n[0]].seq)
			scaff = n[0]
		i += 1
		if sequence:
			outfile2.write('chrUn_random\t'+str(debut+1)+'\t'+str(len(sequence))+'\t'+str(i)+'\tW\t'+scaff+'\t1\t'+str(len(sequence)-debut)+'\t+\n')
			SeqIO.write(SeqRecord(Seq(sequence, generic_dna), id = 'chrUn_random', description=''),outfile, "fasta")
	outfile.close()
	outfile2.close()

if __name__ == "__main__": __main__()