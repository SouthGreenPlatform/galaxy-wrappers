#!/usr/local/bioinfo/python/2.7.9/bin/python
#
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
	for line in file:
		data = line.split()
		if data:
			if data[0][0] == '>':
				if data[0] in dico:
					sys.exit('Two new scaffold '+data[0]+' have the same name')
				else:
					nom = data[0]
					dico[data[0]] = set()
			else:
				for n in dico:
					if data[0] in dico[n]:
						sys.exit('The scaffold '+data[0]+' is already used')
				dico[nom].add(data[0])
	file.close()
	for n in dico:
		if not(n[1:] in dico[n]):
			sys.exit('The new scaffold name '+n[1:]+' is not a used scaffold')

def scaff(TABLE, SEQ, OUT, OUT_VERIF):
	record_dict = SeqIO.index(SEQ, "fasta")
	file = open(TABLE)
	outfile = open(OUT,'w')
	outfile2 = open(OUT_VERIF,'w')
	dico_fait = set()
	sequence = ''
	for line in file:
		data = line.split()
		if data:
			if data[0][0] == '>':
				if sequence:
					SeqIO.write(SeqRecord(Seq(sequence, generic_dna), id = nom, description=''),outfile, "fasta")
					outfile2.write(mot+'\n')
				nom = data[0][1:]
				mot = nom+'\t'
				debut = 1
				sequence = ''
			else:
				if not(data[0] in record_dict):
					sys.exit('The scaffold '+data[0]+' is not in the multifasta')
				if debut:
					debut = 0
					mot = mot+data[0]+'\t'+str(len(sequence)+1)
					if data[1] == "FWD":
						sequence = sequence + str(record_dict[data[0]].seq)
					elif data[1] == "REV":
						sequence = sequence + rev_seq(str(record_dict[data[0]].seq))
					else:
						sys.exit('Orientation information is missing')
					mot = mot+'\t'+str(len(sequence))
				else:
					sequence = sequence + 'NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN'
					mot = mot+'\t'+data[0]+'\t'+str(len(sequence)+1)
					if data[1] == "FWD":
						sequence = sequence + str(record_dict[data[0]].seq)
					elif data[1] == "REV":
						sequence = sequence + rev_seq(str(record_dict[data[0]].seq))
					else:
						sys.exit('Orientation information is missing')
					mot = mot+'\t'+str(len(sequence))
				dico_fait.add(data[0])
	if sequence:
		SeqIO.write(SeqRecord(Seq(sequence, generic_dna), id = nom, description=''),outfile, "fasta")
		outfile2.write(mot+'\n')
	else:
		sys.exit('No sequence in the last scaffold')
	outfile.close()
	return dico_fait

def __main__():
	#Parse Command Line
	parser = optparse.OptionParser()
	parser = optparse.OptionParser(usage="python %prog [options]\n\nProgram designed by Guillaume MARTIN : guillaume.martin@cirad.fr\n\n"
	"This script create junctions between scaffolds using a tabulated file.\n"
	"The input tabulated file look as followed:\n"
	">chr1\n"
	"scaffold1	FWD\n"
	"scaffold2	FWD\n"
	"scaffold3	REV\n"
	">...\n")
	# Wrapper options. 
	parser.add_option( '', '--table', dest='table', default='not_filled', help='The table file of scaffold to join')
	parser.add_option( '', '--fasta', dest='fasta', default='not_filled', help='The multi-fasta scaffold file')
	parser.add_option( '', '--out', dest='out', default='super_contig.fasta', help='The multi-fasta output file name, [default: %default]')
	parser.add_option( '', '--out_verif', dest='out_verif', default='contig2verif.txt', help='The output file to give to verif_fusion.py, [default: %default]')
	(options, args) = parser.parse_args()
	
	
	
	#verifying file	
	verif(options.table)
	
	#creating the scaffolds
	dico_fait = scaff(options.table, options.fasta, options.out, options.out_verif)
	
	#printing the remaining scaffold
	record_dict = SeqIO.index(options.fasta, "fasta")
	outfile = open(options.out,'a')
	for n in record_dict:
		if not(n in dico_fait):
			SeqIO.write(record_dict[n], outfile, "fasta")
	outfile.close()

if __name__ == "__main__": __main__()