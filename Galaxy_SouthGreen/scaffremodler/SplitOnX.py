
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

def __main__():
	#Parse Command Line
	parser = optparse.OptionParser(usage="python %prog [options]\n\nProgram designed by Guillaume MARTIN : guillaume.martin@cirad.fr\n\n"
	"This program split DNA sequence when X are found. All scaffolds are renamed by length.")
	# Wrapper options.
	parser.add_option( '', '--fasta', dest='fasta', default='not_filled', help='The multifasta sequence file')
	parser.add_option( '', '--out', dest='out', default='Splitted.fasta', help='The output file name, [default: %default]')
	(options, args) = parser.parse_args()
	
	
	if options.fasta == 'not_filled':
		sys.exit('--fasta argument is missing')
	
	#loading sequences
	record_dict = SeqIO.index(options.fasta, "fasta")
	
	dico = {}
	liste_taille = []
	i = 0
	for n in record_dict:
		sequence = str(record_dict[n].seq).replace('X',' ').split()
		if len(sequence) == 1:
			i += 1
			if i in dico:
				sys.exit('there is a bug')
			dico[i] = sequence[0]
			liste_taille.append([i, len(sequence[0])])
		else:
			print 'The sequence', n, 'has been cuted', len(sequence)-1, 'time'
			for k in sequence:
				i += 1
				if i in dico:
					sys.exit('there is a bug')
				dico[i] = k
				liste_taille.append([i, len(k)])
				if k[0] == 'N' or k[0] == 'n' or k[0] == 'X':
					sys.exit('Problem at the begining of the sequence : N are found')
				if k[-1] == 'N' or k[-1] == 'n' or k[-1] == 'X':
					sys.exit('Problem at the end of the sequence: N are found')
	
	liste_sorted = sorted(liste_taille, key=operator.itemgetter(1), reverse=True)
	
	outfile = open(options.out,'w')
	i = 0
	for n in liste_sorted:
		i += 1
		SeqIO.write(SeqRecord(Seq(dico[n[0]], generic_dna), id = 'scaffold'+str(i), description=''),outfile, "fasta")
	outfile.close()

if __name__ == "__main__": __main__()