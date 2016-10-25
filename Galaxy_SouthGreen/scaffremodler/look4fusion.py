
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

def stop_err( msg ):
    sys.stderr.write( "%s\n" % msg )
    sys.exit()

def run_job (cmd_line, ERROR):
	print cmd_line
	try:
		tmp = (tempfile.NamedTemporaryFile().name)+'.error'
		# print tmp
		error = open(tmp, 'w')
		proc = subprocess.Popen( args=cmd_line, shell=True, stderr=error)
		returncode = proc.wait()
		error.close()
		error = open( tmp, 'rb' )
		stderr = ''
		buffsize = 1048576
		try:
			while True:
				stderr += error.read( buffsize )
				if not stderr or len( stderr ) % buffsize != 0:
					break
		except OverflowError:
			pass
		error.close()
		os.remove(tmp)
		if returncode != 0:
			raise Exception, stderr
	except Exception, e:
		stop_err( ERROR + str( e ) )


def cherche_fusion(FILE, ZONE, FICHIER, OUT, BOUND):
	outfile = open(OUT, 'a')
	file = open(FILE)
	dico_deb = []
	dico_fin = []
	for line in file:
		data = line.split()
		if data != []:
			# print data
			if int(data[2]) <= ZONE:
				dico_deb.append(data)
			if int(data[1]) >= int(data[3]) - ZONE:
				dico_fin.append(data)
	for n in dico_deb:#on cherche les fusions
		for k in dico_fin:
			if n[0] != k[0]:
				sys.exit('There is a bug')
			if n[6] == k[6]:
				if n[14] == 'AP' and k[14] == 'AP': #destination is before the scaffold searched
					if k[13] == 'BLUE' and n[13] == 'RED':
						if int(n[7])-100 < int(k[7])+100:
							liste_fusion = cherche_N_zone(FICHIER, n[6], n[8], k[7])
							if len(liste_fusion) == 0:
								outfile.write('\t'.join([str(n[0]), str(1), str(n[3]), str(n[4]), str(n[6]), str(n[8]), str(k[7]), str(n[10]), str(n[12]), 'removed : no N found'])+'\n')
							# elif len(liste_fusion) > 1:
								# print 'removed : too much positions', n[0], 0, n[3], n[4], n[6], n[8], k[7], n[10], n[12]
							else:
								for w in liste_fusion:
									outfile.write('\t'.join([str(n[0]), str(1), str(n[3]), str(n[4]), str(n[6]), str(w[0]), str(w[1]), str(n[10]), str(n[12]), 'fusion'])+'\n')
					if k[13] == 'PURPLE' and n[13] == 'GREEN':
						if int(k[7])-100 < int(n[7])+100:
							liste_fusion = cherche_N_zone(FICHIER, n[6], k[8], n[7])
							if len(liste_fusion) == 0:
								outfile.write('\t'.join([str(n[0]), str(1), str(n[3]), str(n[4]), str(n[6]), str(k[8]), str(n[7]), str(n[10]), str(n[12]), 'removed : no N found'])+'\n')
							# elif len(liste_fusion) > 1:
								# print 'removed : too much positions', n[0], 0, n[3], n[4], n[6], k[8], n[7], n[10], n[12]
							else:
								for w in liste_fusion:
									outfile.write('\t'.join([str(n[0]), str(1), str(n[3]), str(n[4]), str(n[6]), str(w[0]), str(w[1]), str(n[10]), str(n[12]), 'fusion'])+'\n')
				elif n[14] == 'AV' and k[14] == 'AV': #destination is before the scaffold searched
					if n[13] == 'BLUE' and k[13] == 'RED':
						if int(n[7])-100 < int(k[7])+100:
							liste_fusion = cherche_N_zone(FICHIER, n[6], n[8], k[7])
							if len(liste_fusion) == 0:
								outfile.write('\t'.join([str(n[0]), str(1), str(n[3]), str(n[4]), str(n[6]), str(n[8]), str(k[7]), str(n[10]), str(n[12]), 'removed : no N found'])+'\n')
							# elif len(liste_fusion) > 1:
								# print 'removed : too much positions', n[0], 0, n[3], n[4], n[6], n[8], k[7], n[10], n[12]
							else:
								for w in liste_fusion:
									outfile.write('\t'.join([str(n[0]), str(1), str(n[3]), str(n[4]), str(n[6]), str(w[0]), str(w[1]), str(n[10]), str(n[12]), 'fusion'])+'\n')
					if k[13] == 'PURPLE' and n[13] == 'GREEN':
						if int(k[7])-100 < int(n[7])+100:
							liste_fusion = cherche_N_zone(FICHIER, n[6], k[8], n[7])
							if len(liste_fusion) == 0:
								outfile.write('\t'.join([str(n[0]), str(1), str(n[3]), str(n[4]), str(n[6]), str(k[8]), str(n[7]), str(n[10]), str(n[12]), 'removed : no N found'])+'\n')
							# elif len(liste_fusion) > 1:
								# print 'removed : too much positions', n[0], 0, n[3], n[4], n[6], k[8], n[7], n[10], n[12]
							else:
								for w in liste_fusion:
									outfile.write('\t'.join([str(n[0]), str(1), str(n[3]), str(n[4]), str(n[6]), str(w[0]), str(w[1]), str(n[10]), str(n[12]), 'fusion'])+'\n')
	#on cherche les contigs
	for n in dico_deb:
		if n[14] == 'AP':
			if n[13] == 'RED':
				if int(n[7]) >= int(n[9]) - ZONE:
					outfile.write('\t'.join([str(n[0]), str(1), str(n[3]), str(n[4]), str(n[6]), str(n[9]), str(n[9]), str(n[10]), str(n[12]), 'contig'])+'\n')
			elif n[13] == 'GREEN':
				if int(n[8]) <= ZONE:
					outfile.write('\t'.join([str(n[0]), str(1), str(n[3]), str(n[4]), str(n[6]), str(1), str(1), str(n[10]), str(n[12]), 'contig'])+'\n')
	for n in dico_fin:
		if n[14] == 'AP':
			if n[13] == 'PURPLE':
				if int(n[7]) >= int(n[9]) - ZONE:
					outfile.write('\t'.join([str(n[0]), str(1), str(n[3]), str(n[4]), str(n[6]), str(n[9]), str(n[9]), str(n[10]), str(n[12]), 'contig'])+'\n')
			elif n[13] == 'BLUE':
				if int(n[8]) <= ZONE:
					outfile.write('\t'.join([str(n[0]), str(1), str(n[3]), str(n[4]), str(n[6]), str(1), str(1), str(n[10]), str(n[12]), 'contig'])+'\n')
	#on cherche des fusion probables
	for n in dico_deb:
		if n[14] == 'AP':
			if n[13] == 'RED':
				liste_fusion2 = cherche_N_zone2(FICHIER, n[6], n[8], 'apres', BOUND)
				if len(liste_fusion2) > 0:
					for w in liste_fusion2:
						outfile.write('\t'.join([str(n[0]), str(1), str(n[3]), str(n[4]), str(n[6]), str(w[0]), str(w[1]), str(n[10]), str(n[12]), 'prob_fusion'])+'\n')
			if n[13] == 'GREEN':
				liste_fusion2 = cherche_N_zone2(FICHIER, n[6], n[7], 'avant', BOUND)
				if len(liste_fusion2) > 0:
					for w in liste_fusion2:
						outfile.write('\t'.join([str(n[0]), str(1), str(n[3]), str(n[4]), str(n[6]), str(w[0]), str(w[1]), str(n[10]), str(n[12]), 'prob_fusion'])+'\n')
		if n[14] == 'AV':
			if n[13] == 'BLUE':
				liste_fusion2 = cherche_N_zone2(FICHIER, n[6], n[8], 'apres', BOUND)
				if len(liste_fusion2) > 0:
					for w in liste_fusion2:
						outfile.write('\t'.join([str(n[0]), str(1), str(n[3]), str(n[4]), str(n[6]), str(w[0]), str(w[1]), str(n[10]), str(n[12]), 'prob_fusion'])+'\n')
	for n in dico_fin:
		if n[14] == 'AP':
			if n[13] == 'PURPLE':
				liste_fusion2 = cherche_N_zone2(FICHIER, n[6], n[8], 'apres', BOUND)
				if len(liste_fusion2) > 0:
					for w in liste_fusion2:
						outfile.write('\t'.join([str(n[0]), str(1), str(n[3]), str(n[4]), str(n[6]), str(w[0]), str(w[1]), str(n[10]), str(n[12]), 'prob_fusion'])+'\n')
			if n[13] == 'BLUE':
				liste_fusion2 = cherche_N_zone2(FICHIER, n[6], n[7], 'avant', BOUND)
				if len(liste_fusion2) > 0:
					for w in liste_fusion2:
						outfile.write('\t'.join([str(n[0]), str(1), str(n[3]), str(n[4]), str(n[6]), str(w[0]), str(w[1]), str(n[10]), str(n[12]), 'prob_fusion'])+'\n')
		if n[14] == 'AV':
			if n[13] == 'RED':
				liste_fusion2 = cherche_N_zone2(FICHIER, n[6], n[7], 'avant', BOUND)
				if len(liste_fusion2) > 0:
					for w in liste_fusion2:
						outfile.write('\t'.join([str(n[0]), str(1), str(n[3]), str(n[4]), str(n[6]), str(w[0]), str(w[1]), str(n[10]), str(n[12]), 'prob_fusion'])+'\n')
	outfile.close()

def cherche_N_zone2(FILE, SCAFFOLD, POS, WHERE, MAX_DIST):
	liste = []
	file = open(FILE)
	debut = ''
	fin = ''
	for line in file:
		data = line.split()
		if data:
			info = data[3].split('-')
			if len(info) == 3:
				if info[0] == SCAFFOLD:
					if WHERE == 'avant':
						if info[2] == 'v':
							debut = int(info[1])
						elif info[2] == 'p' and debut and (int(POS) + 100) >= int(info[1]) and int(POS) <= (int(info[1]) + MAX_DIST):
							liste.append([debut,int(info[1])])
							debut = ''
					elif WHERE == 'apres':
						if info[2] == 'v' and (int(POS) - 100) <= int(info[1]) and int(POS) >= (int(info[1]) - MAX_DIST):
							debut = int(info[1])
						elif info[2] == 'p' and debut:
							liste.append([debut,int(info[1])])
							debut = ''
			
	return liste
						

def cherche_N_zone(FILE, SCAFFOLD, DEBUT, FIN):
	liste = []
	file = open(FILE)
	debut = ''
	for line in file:
		data = line.split()
		if data:
			info = data[3].split('-')
			if len(info) == 3:
				if info[0] == SCAFFOLD:
					if (int(DEBUT) - 100) <= int(info[1]) and int(info[1]) <= (int(FIN) + 100):
						if info[2] == 'v':
							debut = int(info[1])
						elif info[2] == 'p' and debut:
							liste.append([debut,int(info[1])])
							debut = ''
						else:
							debut = ''
	return liste

def look4fusion(LOCA_PROGRAMS, CONFIG, BOUND, OUT_TEXT, OUT_TAR, PATHNAME):
	config = ConfigParser.RawConfigParser()
	config.read(CONFIG)
	
	dic_chr = {}
	chr_order = []
	file = open(config.get('General','chr'))
	for line in file:
		data = line.split()
		if data:
			chr_order.append(data[0])
			dic_chr[data[0]] = int(data[1])
	file.close()
	
	
	liste_out = []
	file = open(config.get('General','chr'))
	for line in file:
		data = line.split()
		liste_region = []
		if data:
			if data[0][0] != "#":
				tempo = tempfile.NamedTemporaryFile()
				print tempo.name
				if config.get('General','orient') == 'rf':
					if config.get('General','chr_fr') == 'yes':
						if int(data[1]) < BOUND:
							liste_region = liste_region + cherche_dest(data[0], 1, int(data[1]), config.get('Discord_zone','chr_fr'), tempo, 'BLUE', 'FWD', chr_order, dic_chr)
						else:
							liste_region = liste_region + cherche_dest(data[0], 1, BOUND, config.get('Discord_zone','chr_fr'), tempo, 'BLUE', 'FWD', chr_order, dic_chr)
							liste_region = liste_region + cherche_dest(data[0], int(data[1])-BOUND, int(data[1]), config.get('Discord_zone','chr_fr'), tempo, 'BLUE', 'FWD', chr_order, dic_chr)
					if config.get('General','chr_rf') == 'yes':
						if int(data[1]) < BOUND:
							liste_region = liste_region + cherche_dest(data[0], 1, int(data[1]), config.get('Discord_zone','chr_rf'), tempo, 'RED', 'FWD', chr_order, dic_chr)
						else:
							liste_region = liste_region + cherche_dest(data[0], 1, BOUND, config.get('Discord_zone','chr_rf'), tempo, 'RED', 'FWD', chr_order, dic_chr)
							liste_region = liste_region + cherche_dest(data[0], int(data[1])-BOUND, int(data[1]), config.get('Discord_zone','chr_rf'), tempo, 'RED', 'FWD', chr_order, dic_chr)
					if config.get('General','chr_ff') == 'yes':
						if int(data[1]) < BOUND:
							liste_region = liste_region + cherche_dest(data[0], 1, int(data[1]), config.get('Discord_zone','chr_ff'), tempo, 'GREEN', 'REV', chr_order, dic_chr)
						else:
							liste_region = liste_region + cherche_dest(data[0], 1, BOUND, config.get('Discord_zone','chr_ff'), tempo, 'GREEN', 'REV', chr_order, dic_chr)
							liste_region = liste_region + cherche_dest(data[0], int(data[1])-BOUND, int(data[1]), config.get('Discord_zone','chr_ff'), tempo, 'GREEN', 'REV', chr_order, dic_chr)
					if config.get('General','chr_rr') == 'yes':
						if int(data[1]) < BOUND:
							liste_region = liste_region + cherche_dest(data[0], 1, int(data[1]), config.get('Discord_zone','chr_rr'), tempo, 'PURPLE', 'REV', chr_order, dic_chr)
						else:
							liste_region = liste_region + cherche_dest(data[0], 1, BOUND, config.get('Discord_zone','chr_rr'), tempo, 'PURPLE', 'REV', chr_order, dic_chr)
							liste_region = liste_region + cherche_dest(data[0], int(data[1])-BOUND, int(data[1]), config.get('Discord_zone','chr_rr'), tempo, 'PURPLE', 'REV', chr_order, dic_chr)
				elif config.get('General','orient') == 'fr':
					if config.get('General','chr_fr') == 'yes':
						if int(data[1]) < BOUND:
							liste_region = liste_region + cherche_dest(data[0], 1, int(data[1]), config.get('Discord_zone','chr_fr'), tempo, 'RED', 'FWD', chr_order, dic_chr)
						else:
							liste_region = liste_region + cherche_dest(data[0], 1, BOUND, config.get('Discord_zone','chr_fr'), tempo, 'RED', 'FWD', chr_order, dic_chr)
							liste_region = liste_region + cherche_dest(data[0], int(data[1])-BOUND, int(data[1]), config.get('Discord_zone','chr_fr'), tempo, 'RED', 'FWD', chr_order, dic_chr)
					if config.get('General','chr_rf') == 'yes':
						if int(data[1]) < BOUND:
							liste_region = liste_region + cherche_dest(data[0], 1, int(data[1]), config.get('Discord_zone','chr_rf'), tempo, 'BLUE', 'FWD', chr_order, dic_chr)
						else:
							liste_region = liste_region + cherche_dest(data[0], 1, BOUND, config.get('Discord_zone','chr_rf'), tempo, 'BLUE', 'FWD', chr_order, dic_chr)
							liste_region = liste_region + cherche_dest(data[0], int(data[1])-BOUND, int(data[1]), config.get('Discord_zone','chr_rf'), tempo, 'BLUE', 'FWD', chr_order, dic_chr)
					if config.get('General','chr_ff') == 'yes':
						if int(data[1]) < BOUND:
							liste_region = liste_region + cherche_dest(data[0], 1, int(data[1]), config.get('Discord_zone','chr_ff'), tempo, 'PURPLE', 'REV', chr_order, dic_chr)
						else:
							liste_region = liste_region + cherche_dest(data[0], 1, BOUND, config.get('Discord_zone','chr_ff'), tempo, 'PURPLE', 'REV', chr_order, dic_chr)
							liste_region = liste_region + cherche_dest(data[0], int(data[1])-BOUND, int(data[1]), config.get('Discord_zone','chr_ff'), tempo, 'PURPLE', 'REV', chr_order, dic_chr)
					if config.get('General','chr_rr') == 'yes':
						if int(data[1]) < BOUND:
							liste_region = liste_region + cherche_dest(data[0], 1, int(data[1]), config.get('Discord_zone','chr_rr'), tempo, 'GREEN', 'REV', chr_order, dic_chr)
						else:
							liste_region = liste_region + cherche_dest(data[0], 1, BOUND, config.get('Discord_zone','chr_rr'), tempo, 'GREEN', 'REV', chr_order, dic_chr)
							liste_region = liste_region + cherche_dest(data[0], int(data[1])-BOUND, int(data[1]), config.get('Discord_zone','chr_rr'), tempo, 'GREEN', 'REV', chr_order, dic_chr)
				else:
					sys.exit('Unsupported orientation found in the conf file')
				tempo.flush()
				cherche_fusion(tempo.name, BOUND, config.get('General','out_n'), OUT_TEXT, BOUND)
				tempo.close()
				if liste_region:
					if (BOUND*2) >= int(data[1]):
						mot = data[0]+':0:'+str(int(data[1]))
					else:
						mot = data[0]+':0:'+str(BOUND)
						mot = mot+'-'+data[0]+':'+str(int(data[1])-int(BOUND))+':'+str(int(data[1]))
					for n in liste_region:
						if ((n[1] - BOUND) <= 0) and ((n[2] + BOUND) >= dic_chr[n[0]]):
							mot = mot+'-'+n[0]+':0:'+str(dic_chr[n[0]])
						elif ((n[1] - BOUND) <= 0):
							mot = mot+'-'+n[0]+':0:'+str(n[2] + BOUND)
						elif((n[2] + BOUND) >= dic_chr[n[0]]):
							mot = mot+'-'+n[0]+':'+str(n[1] - BOUND)+':'+str(dic_chr[n[0]])
						else:
							mot = mot+'-'+n[0]+':'+str(n[1] - BOUND)+':'+str(n[2] + BOUND)
					liste_out.append(data[0]+'.png')
					redraw = '%s %s/draw_circos.py --config %s --frf n --ff n --rr n --ins n --delet n --chr_rr n --chr_rf n --chr_fr n --chr_ff n --draw %s --out %s' % (LOCA_PROGRAMS.get('Programs','python'), PATHNAME, CONFIG, mot, data[0]+'.png')
					run_job(redraw, 'Bug when drawing circos')
	mot = liste_out[0]
	for n in liste_out[1:]:
		mot = mot +' '+ n
	archivage = 'tar -cf '+OUT_TAR+' '+mot
	run_job(archivage, 'Bug in archive creation')
	for n in liste_out:
		os.remove(n)

def cherche_dest(CHR, DEBUT, FIN, FICHIER, TEMP, COLOR, TYPE, CHR_ORDER, DIC_CHR):
	liste = []
	file = open(FICHIER)
	for line in file:
		data = line.split()
		if data:
			# print data
			if data[0] != data[5]  and data[13] == "PASSED":
				if CHR == data[0]:
					if DEBUT-100 <= int(data[1]) and int(data[2]) <= FIN+100:
						liste.append([data[5], int(data[6]), int(data[7])])
						if CHR_ORDER.index(data[0]) < CHR_ORDER.index(data[5]):
							POSITION = 'AV'
						elif CHR_ORDER.index(data[0]) > CHR_ORDER.index(data[5]):
							POSITION = 'AP'
						else:
							sys.exit('bug in cherche_dest')
						TEMP.write('\t'.join([data[0], data[1], data[2], str(DIC_CHR[data[0]]), 'FWD', '-', data[5], data[6], data[7], str(DIC_CHR[data[5]]), 'FWD', '-', TYPE, COLOR, POSITION])+'\n')
				elif CHR == data[5]:
					if DEBUT-100 <= int(data[6]) and int(data[7]) <= FIN+100:
						liste.append([data[0], int(data[1]), int(data[2])])
						if CHR_ORDER.index(data[5]) < CHR_ORDER.index(data[0]):
							POSITION = 'AV'
						elif CHR_ORDER.index(data[5]) > CHR_ORDER.index(data[0]):
							POSITION = 'AP'
						else:
							sys.exit('bug in cherche_dest')
						TEMP.write('\t'.join([data[5], data[6], data[7], str(DIC_CHR[data[5]]), 'FWD', '-', data[0], data[1], data[2], str(DIC_CHR[data[0]]), 'FWD', '-', TYPE, COLOR, POSITION])+'\n')
	TEMP.flush()
	return liste

def __main__():
	#Parse Command Line
	parser = optparse.OptionParser(usage="python %prog [options]\n\nProgram designed by Guillaume MARTIN : guillaume.martin@cirad.fr\n\n"
	"This script looks for possible scaffold fusions and junctions")
	# Wrapper options. 
	parser.add_option( '', '--config', dest='config', default='not_filled', help='The conf file generated by conf4circos.py')
	
	parser.add_option( '', '--bound', dest='bound', default=10000, help='Boudary of scaffold to look for fusion and junction, [default: %default]')
	parser.add_option( '', '--out', dest='out', default='possible_fusion.txt', help='Output text file of possible fusion and contig, [default: %default]')
	parser.add_option( '', '--out_tar', dest='out_tar', default='possible_fusion.tar', help='Output tar.gz file containing circos figures, [default: %default]')
	
	
	(options, args) = parser.parse_args()
	
	pathname = os.path.dirname(sys.argv[0])
	
	loca_programs = ConfigParser.RawConfigParser()
	loca_programs.read(pathname+'/loca_programs.conf')
	
	temp = tempfile.NamedTemporaryFile().name
	print temp
	
	config = ConfigParser.RawConfigParser()
	config.read(options.config)
	
	if options.config == 'not_filled':
		sys.exit('--config argument is missing')
	if config.get('General','orient') != 'rf':
		sys.exit('The program exited: only rf orientation is accepted')
	
	look4fusion(loca_programs, options.config, int(options.bound), temp, options.out_tar, pathname)
	
	os.system('cat '+temp+' | sort | uniq > '+options.out)
	os.remove(temp)

if __name__ == "__main__": __main__()