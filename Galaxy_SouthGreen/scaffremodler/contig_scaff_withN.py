
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

def verif(TABLE):
	dico = {}
	file = open(TABLE)
	prec = ''
	for line in file:
		data = line.split()
		if data:
			if data[0] == prec:
				if data[1] == 'W':
					for n in dico:
						if data[2] in dico[n]:
							mot = 'Warning: the scaffold '+data[2]+' is already used'
							print mot
						dico[nom].add(data[0])
			else:
				if data[0] in dico:
					mot = 'Two new scaffold '+data[0]+' have the same name'
					sys.exit(mot)
				else:
					nom = data[0]
					dico[data[0]] = set()
				prec = data[0]
	file.close()
	for n in dico:
		if not(n in dico[n]):
			sys.exit('The new scaffold name '+n[1:]+' is not a used scaffold')

def scaff(TABLE, SEQ, OUT, OUT_VERIF):
	record_dict = SeqIO.index(SEQ, "fasta")
	file = open(TABLE)
	outfile = open(OUT,'w')
	outfile2 = open(OUT_VERIF,'w')
	dico_fait = set()
	sequence = ''
	prec = ''
	for line in file:
		data = line.split()
		if data:
			if data[0] != prec:
				if sequence:
					SeqIO.write(SeqRecord(Seq(sequence, generic_dna), id = prec, description=''),outfile, "fasta")
				prec = data[0]
				mot = data[0]+'\t'
				j = 0
				sequence = ''
				if data[1] == 'W':
					debut = len(sequence)
					j += 1
					if data[3] == '+':
						sequence = sequence + str(record_dict[data[2]].seq)
					elif data[3] == '-':
						sequence = sequence + rev_seq(str(record_dict[data[2]].seq))
					else:
						sys.exit('Wrong orientation information'+data[3])
					outfile2.write(prec+'\t'+str(debut+1)+'\t'+str(len(sequence))+'\t'+str(j)+'\tW\t'+data[2]+'\t1\t'+str(len(sequence)-debut)+'\t+\n')
					dico_fait.add(data[2])
				elif data[1] == 'N':
					j += 1
					i = 0
					N_number = int(data[2])
					debut = len(sequence)
					while i < N_number:
						sequence = sequence + 'N'
						i += 1
					outfile2.write(prec+'\t'+str(debut+1)+'\t'+str(len(sequence))+'\t'+str(j)+'\tN\t'+data[2]+'\tfragment\tno\n')
				else:
					sys.exit('Wrong region type'+data[3])
			else:
				if data[1] == 'W':
					debut = len(sequence)
					j += 1
					if data[3] == '+':
						sequence = sequence + str(record_dict[data[2]].seq)
					elif data[3] == '-':
						sequence = sequence + rev_seq(str(record_dict[data[2]].seq))
					else:
						sys.exit('Wrong orientation information'+data[3])
					outfile2.write(prec+'\t'+str(debut+1)+'\t'+str(len(sequence))+'\t'+str(j)+'\tW\t'+data[2]+'\t1\t'+str(len(sequence)-debut)+'\t+\n')
					dico_fait.add(data[2])
				elif data[1] == 'N':
					j += 1
					i = 0
					N_number = int(data[2])
					debut = len(sequence)
					while i < N_number:
						sequence = sequence + 'N'
						i += 1
					outfile2.write(prec+'\t'+str(debut+1)+'\t'+str(len(sequence))+'\t'+str(j)+'\tN\t'+data[2]+'\tfragment\tno\n')
				else:
					sys.exit('Wrong region type'+data[3])
	if sequence:
		SeqIO.write(SeqRecord(Seq(sequence, generic_dna), id = prec, description=''),outfile, "fasta")
	else:
		sys.exit('No sequence in the last scaffold')
	outfile.close()
	return dico_fait

def __main__():
	#Parse Command Line
	parser = optparse.OptionParser()
	parser = optparse.OptionParser(usage="python %prog [options]\n\nProgram designed by Guillaume MARTIN : guillaume.martin@cirad.fr\n\n"
	"This script create junctions between scaffolds using a tabulated file.\n")
	# Wrapper options. 
	parser.add_option( '', '--table', dest='table', default='not_filled', help='The table file of scaffold to join')
	parser.add_option( '', '--fasta', dest='fasta', default='not_filled', help='The multi-fasta scaffold file')
	parser.add_option( '', '--out', dest='out', default='super_contig.fasta', help='The multi-fasta output file name, [default: %default]')
	parser.add_option( '', '--out_info', dest='out_info', default='contig_info.agp', help='An agp file locating contigs in scaffold, [default: %default]')
	(options, args) = parser.parse_args()
	
	
	
	
	#verifying file	
	verif(options.table)
	
	#creating the scaffolds
	dico_fait = scaff(options.table, options.fasta, options.out, options.out_info)
	
	#printing the remaining scaffold
	record_dict = SeqIO.index(options.fasta, "fasta")
	outfile = open(options.out,'a')
	for n in record_dict:
		if not(n in dico_fait):
			SeqIO.write(record_dict[n], outfile, "fasta")
	outfile.close()

if __name__ == "__main__": __main__()