
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
		dic[liste[0]] = liste[2]
	return dic


def Filtre(SAM, ASXS, QUAL, OUT):
	outfile = open(OUT, 'w')
	
	nb_input = 0
	nb_kept = 0
	
	file = open(SAM)
	if QUAL == 'not_filled':
		min_dif = int(ASXS)
		l1 = file.readline().split()
		while l1[0][0] == '@':
			outfile.write('\t'.join(l1)+'\n')
			l1 = file.readline().split()
		l2 = file.readline().split()
		while l1:
			nb_input += 1
			if l1[0] != l2[0]:
				sys.exit('Read should be sorted by query name in the sam file')
			dico1 = find_info(l1)
			dico2 = find_info(l2)
			if 'XS' in dico1:
				if 'XS' in dico2:
					if abs(int(dico1['AS'])-int(dico1['XS'])) >= min_dif and abs(int(dico2['AS'])-int(dico2['XS'])) >= min_dif:
						outfile.write('\t'.join(l1)+'\n')
						outfile.write('\t'.join(l2)+'\n')
						nb_kept += 1
				else:
					if abs(int(dico1['AS'])-int(dico1['XS'])) >= min_dif:
						outfile.write('\t'.join(l1)+'\n')
						outfile.write('\t'.join(l2)+'\n')
						nb_kept += 1
			elif 'XS' in dico2:
				if abs(int(dico2['AS'])-int(dico2['XS'])) >= min_dif:
					outfile.write('\t'.join(l1)+'\n')
					outfile.write('\t'.join(l2)+'\n')
					nb_kept += 1
			else:
				outfile.write('\t'.join(l1)+'\n')
				outfile.write('\t'.join(l2)+'\n')
				nb_kept += 1
			l1 = file.readline().split()
			l2 = file.readline().split()

	elif ASXS == 'not_filled':
		min_qual = int(QUAL)
		l1 = file.readline().split()
		while l1[0][0] == '@':
			outfile.write('\t'.join(l1)+'\n')
			l1 = file.readline().split()
		l2 = file.readline().split()
		while l1:
			nb_input += 1
			if l1[0] != l2[0]:
				sys.exit('Read should be sorted by query name in the sam file')
			if int(l1[4]) >= min_qual and int(l2[4]) >= min_qual:
				outfile.write('\t'.join(l1)+'\n')
				outfile.write('\t'.join(l2)+'\n')
				nb_kept += 1
			l1 = file.readline().split()
			l2 = file.readline().split()
	else:
		min_dif = int(ASXS)
		min_qual = int(QUAL)
		l1 = file.readline().split()
		while l1[0][0] == '@':
			outfile.write('\t'.join(l1)+'\n')
			l1 = file.readline().split()
		l2 = file.readline().split()
		while l1:
			nb_input += 1
			if l1[0] != l2[0]:
				sys.exit('Read should be sorted by query name in the sam file')
			dico1 = find_info(l1)
			dico2 = find_info(l2)
			if 'XS' in dico1:
				if 'XS' in dico2:
					if abs(int(dico1['AS'])-int(dico1['XS'])) >= min_dif and abs(int(dico2['AS'])-int(dico2['XS'])) >= min_dif:
						if int(l1[4]) >= min_qual and int(l2[4]) >= min_qual:
							outfile.write('\t'.join(l1)+'\n')
							outfile.write('\t'.join(l2)+'\n')
							nb_kept += 1
				else:
					if abs(int(dico1['AS'])-int(dico1['XS'])) >= min_dif:
						if int(l1[4]) >= min_qual and int(l2[4]) >= min_qual:
							outfile.write('\t'.join(l1)+'\n')
							outfile.write('\t'.join(l2)+'\n')
							nb_kept += 1
			elif 'XS' in dico2:
				if abs(int(dico2['AS'])-int(dico2['XS'])) >= min_dif:
					if int(l1[4]) >= min_qual and int(l2[4]) >= min_qual:
						outfile.write('\t'.join(l1)+'\n')
						outfile.write('\t'.join(l2)+'\n')
						nb_kept += 1
			else:
				if int(l1[4]) >= min_qual and int(l2[4]) >= min_qual:
					outfile.write('\t'.join(l1)+'\n')
					outfile.write('\t'.join(l2)+'\n')
					nb_kept += 1
			l1 = file.readline().split()
			l2 = file.readline().split()

	print('Mapped pairs: %s' % nb_input)
	print('Mapped pairs kept: %s' % nb_kept)
	print('Mapped pairs proportion kept: %s' % str(float(nb_kept)/float(nb_input)))

def __main__():
	#Parse Command Line
	parser = optparse.OptionParser(usage="python %prog [options]\n\nProgram designed by Guillaume MARTIN : guillaume.martin@cirad.fr"
	"\n\n This script takes a sam file and output only paired reads were both mates pass filter threshold on the AS/XS flags or mapping quality or both provided.")
	# Wrapper options. 
	parser.add_option( '', '--sam', dest='sam', default='not_filled', help='Paired sam file')
	parser.add_option( '', '--asxs', dest='asxs', default='not_filled', help='Minimal difference between the best and second hit accepted to consider the hit as single')
	parser.add_option( '', '--qual', dest='qual', default='not_filled', help='Minimal mapping quality to keep the hit')
	parser.add_option( '', '--rminput', dest='rminput', default='n', help='Remove input file: y or n, [default: %default]')
	parser.add_option( '', '--out', dest='out', default='Single_hit_mapped.sam', help='Output file')
	parser.add_option( '', '--config', dest='config', default=None)
	(options, args) = parser.parse_args()
	
	
	
	if options.config:
		config = ConfigParser.RawConfigParser()
		config.read(options.config)
		if config.get('Single_filter','filter_multi') == 'y':
			Filtre(config.get('Mapping','out'), config.get('Single_filter','asxs'), config.get('Single_filter','qual'), options.out)
		else:
			print 'The input sam is the sam as the output sam in 3_filter_single_pair'
			os.system('cp % %' % (config.get('Mapping','out'), options.out))
		if config.get('Single_filter','rminput') == 'y':
			os.remove(config.get('Mapping','out'))
		config.set('Single_filter', 'out', options.out)
		config.set('Single_filter', 'type', 'sam')
		config.set('Remove_dup', 'sort', 'coordinate')
		with open(options.config, 'wb') as configfile:
			config.write(configfile)
	else:
		if options.sam == 'not_filled':
			sys.exit('--sam argument is missing')
		if options.qual == 'not_filled' and options.asxs == 'not_filled':
			print 'No --asxs or --qual argument are passed the imput sam is the sam as the output sam in 3_filter_single_pair'
			os.system('cp % %' % (options.sam, options.out))
		else:
			Filtre(options.sam, options.asxs, options.qual, options.out)
		if options.rminput == 'y':
			os.remove(options.sam)


if __name__ == "__main__": __main__()



