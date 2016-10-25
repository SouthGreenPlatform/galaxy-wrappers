
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
	"This program replace specified regions in the provided table file by X. These X will be used to split scaffold using SplitOnX.py"
	"The table file should be formated has in the example:"
	"scaffold83 93565 93586"
	"scaffold120 330181 330183"
	"scaffold120 380870 383428")
	# Wrapper options. 
	parser.add_option( '', '--table', dest='table', default='not_filled', help='The table file with region to convert to X')
	parser.add_option( '', '--fasta', dest='fasta', default='not_filled', help='The multifasta sequence file')
	parser.add_option( '', '--out', dest='out', default='X_converted.fasta', help='The output file name, [default: %default]')
	(options, args) = parser.parse_args()
	
	
	
	if options.table == 'not_filled':
		sys.exit('--table argument is missing')
	if options.fasta == 'not_filled':
		sys.exit('--fasta argument is missing')
	
	#loading sequences
	record_dict = SeqIO.index(options.fasta, "fasta")
	file = open(options.table)
	dic = {}
	for line in file:
		data = line.split()
		if data:
			if data[0] in dic:
				if len(data) == 2:
					dic[data[0]].add(int(data[1])-1)
				else:
					i = int(data[1])
					while i <= int(data[2]):
						dic[data[0]].add(i-1)
						i += 1
			else:
				dic[data[0]] = set()
				if len(data) == 2:
					dic[data[0]].add(int(data[1])-1)
				else:
					i = int(data[1])
					while i <= int(data[2]):
						dic[data[0]].add(i-1)
						i += 1
	file.close()
	
	outfile = open(options.out,'w')
	for n in record_dict:
		if n in dic:
			sequence = list(str(record_dict[n].seq))
			for k in dic[n]:
				sequence[k] = 'X'
			SeqIO.write(SeqRecord(Seq(''.join(sequence), generic_dna), id = n, description=''),outfile, "fasta")
		else:
			SeqIO.write(SeqRecord(record_dict[n].seq, id = n, description=''),outfile, "fasta")

if __name__ == "__main__": __main__()