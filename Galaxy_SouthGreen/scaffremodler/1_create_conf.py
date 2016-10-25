
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

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def cree_chrom(FILE, OUT):
	record_dict = SeqIO.index(FILE, "fasta")
	outfile = open(OUT, 'wb')
	liste = []
	for n in record_dict:
		liste.append(n)
	liste.sort()
	for n in liste:
		outfile.write('\t'.join([n, str(len(str(record_dict[n].seq)))])+'\n')
	outfile.close()

def __main__():
	#Parse Command Line
	parser = optparse.OptionParser(usage="python %prog [options]\n\nProgram designed by Guillaume MARTIN : guillaume.martin@cirad.fr\n\n"
	"This script generate a configuration file that will be used in the ApMap pipeline")
	# Wrapper options.
	parser.add_option( '', '--tool', dest='tool', default='bowtie2_single', help='The tool used : bowtie, bowtie2, bowtie2_single, bwa, bwa_mem, [default: %default]')
	parser.add_option( '', '--ref', dest='ref', default='not_filled', help='The multifasta reference file')
	parser.add_option( '', '--q1', dest='q1', default='not_filled', help='The mate1 fastq file')
	parser.add_option( '', '--q2', dest='q2', default='not_filled', help='The mate2 fastq file')
	parser.add_option( '', '--orient', dest='orient', default='rf', help='The expected orientation: rf or fr, [default: %default]')
	parser.add_option( '', '--mini', dest='mini', default='2500', help='The minimum insert size (integer), [default: %default]')
	parser.add_option( '', '--maxi', dest='maxi', default='7500', help='The maximum insert size (integer), [default: %default]')
	parser.add_option( '', '--qual', dest='qual', default='33', help='Fastq quality encoding: 33 or 64, [default: %default]')
	parser.add_option( '', '--index', dest='index', default='y', help='Build reference index : y or n,  [default: %default]')
	parser.add_option( '', '--rmindex', dest='rmindex', default='y', help='Remove reference index at the end of calculation: y or n, [default: %default]')
	parser.add_option( '', '--filter_multi', dest='filter_multi', default='y', help='Filter reads with multiple locations : y or n,  [default: %default]')
	parser.add_option( '', '--mini_dis', dest='mini_dis', default='10000', help='The minimal insert size to keep the discordant read for structural variation search (integer), [default: %default]')
	parser.add_option( '', '--mult_max_cov', dest='mult_max_cov', default='10', help='multiplicator of median coverage for maximal median coverage to keep a zone (float), [default: %default]')
	parser.add_option( '', '--mult_min_cov', dest='mult_min_cov', default='0.25', help='multiplicator of median coverage for minimal median coverage to keep a zone (float),  [default: %default]')
	parser.add_option( '', '--min_zone', dest='min_zone', default='500', help='Minimal number of covered sites in a zone to be considered (integer),  [default: %default]')
	parser.add_option( '', '--min_gap', dest='min_gap', default='300', help='Maximal number of contiguous uncovered sites in a zone to be considered as a single zone (integer),  [default: %default]')
	parser.add_option( '', '--thread', dest='thread', default='1', help='The thread number used for mapping (integer), [default: %default]')
	parser.add_option( '', '--msd', dest='msd', default='3', help='Multiplicator of standard deviation for discordant zone identification (integer), [default: %default]')
	parser.add_option( '', '--max_dist_merge', dest='max_dist_merge', default=1000, help='Maximal distance between two discordant zone to merge, [default: %default]')
	parser.add_option( '', '--YiS', dest='YiS', default=0, help='The Y-intercept of the linear function for zone size that will give the first component of product giving the score (integer), [default: %default]')
	parser.add_option( '', '--MiS', dest='MiS', default=0.5, help='Multiplicator of median insert size for calculating minimal zone size for which the first component of product giving the score will be maximal (integer), [default: %default]. Exmple: if 0.5, discordant zone of more than 2500 pb will have a maximal score')
	parser.add_option( '', '--YiC', dest='YiC', default=0, help='The Y-intercept of the linear function for coverage that will give the second component of product giving the score (integer), [default: %default]')
	parser.add_option( '', '--MiC', dest='MiC', default=0.25, help='Multiplicator of median coverage for calculating  minimal zone coverage for which the second component of product giving the score will be maximal (integer), [default: %default].  For homozygous SV in diploid: expected value = 0.5, if heterozygous: expected value = 0.25')
	parser.add_option( '', '--min_score', dest='min_score', default=70, help='The minimal score for a discordant zone to be identified as passed, [default: %default]')
	parser.add_option( '', '--ploid', dest='ploid', default=0.33, help='Multiplicator for coverage variation detection in SV identification (ex : If homozygous duplication expected in diploid: expected = coverage + coverage*1, if heterozygous duplication expected in diploid => expected = coverage + coverage*0.5). Choose a value lower than the expected one')
	parser.add_option( '', '--restimate', dest='restimate', default='n', help='Wether re-estimating --mini and --maxi parameters: y or n, [default: %default]. If y, these parameters are calculated as followed on well mapped paired read on the basis of previous min and max parameters: min/max = median -/+ (standard_deviation * "--msd" option)')
	parser.add_option( '', '--output', dest='output', default='config.conf', help='The output of the conf file, [default: %default]')
	parser.add_option( '', '--chr', dest='chr', default='chr.tab', help='Output file containing chromosomes informations, [default: %default]')
	parser.add_option( '', '--rm_intermediate', dest='rm_intermediate', default='n', help='remove intermediate bam/sam, [default: %default]')
	parser.add_option( '', '--exclude_chrom', dest='exclude_chrom', default='no_exclude', help='Exclude chromosomes from analysis. "no_exclude" or chromosomes names separated by "=", [default: %default]')
	(options, args) = parser.parse_args()
	
	
	
	cree_chrom(options.ref, options.chr)
	# print options.ref
	# print options.chr
	# print options.q1
	# print options.q2
	# print options.chr
	# print options.output

	config = ConfigParser.RawConfigParser()
	config.add_section('General')
	config.set('General','ref', options.ref)
	config.set('General','chr', options.chr)
	config.set('General','mini', options.mini)
	config.set('General','maxi', options.maxi)
	config.set('General','thread', options.thread)
	config.set('General','tool', options.tool)
	config.set('General','q1', options.q1)
	config.set('General','q2', options.q2)
	config.set('General','qual', options.qual)
	config.set('General','orient', options.orient)
	config.set('General','index', options.index)
	config.set('General','rmindex', options.rmindex)
	config.set('General','sd_multiplicator', options.msd)
	config.set('General','restimate', options.restimate)
	config.set('General','mini_dis', options.mini_dis)
	config.set('General','mult_max_cov', options.mult_max_cov)
	config.set('General','mult_min_cov', options.mult_min_cov)
	config.set('General','min_zone', options.min_zone)
	config.set('General','min_gap', options.min_gap)
	config.set('General','max_dist_merge', options.max_dist_merge)
	config.set('General','YiS', options.YiS)
	config.set('General','MiS', options.MiS)
	config.set('General','YiC', options.YiC)
	config.set('General','MiC', options.MiC)
	config.set('General','min_score', options.min_score)
	config.set('General','ploid', options.ploid)
	config.set('General','fai_file', options.ref+'.fai')
	config.set('General','exclude_chrom', options.exclude_chrom)
	config.add_section('Mapping')
	config.add_section('Single_filter')
	config.set('Single_filter','rminput', options.rm_intermediate)
	config.set('Single_filter','filter_multi', options.filter_multi)
	config.add_section('Remove_dup')
	config.set('Remove_dup','rminput', options.rm_intermediate)
	config.add_section('Calc_coverage')
	config.add_section('Trie_discord')
	config.set('Trie_discord','rminput', options.rm_intermediate)
	config.add_section('Score_discord')
	config.add_section('Ident_discord')
	with open(options.output, 'wb') as configfile:
		config.write(configfile)
	
if __name__ == "__main__": __main__()
