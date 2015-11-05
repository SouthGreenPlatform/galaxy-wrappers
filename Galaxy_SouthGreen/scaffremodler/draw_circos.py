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

import optparse, os, shutil, subprocess, sys, tempfile, fileinput, ConfigParser, operator, time, math, datetime

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

def define_regions(DRAW, CHR):
	#Record regions to draw
	dic = {}
	chrDic = {}
	taille_total = 0
	f = open(CHR, 'r')
	for line in f:
		if line.strip():
			cols = line.split()
			chrDic[cols[0]] = int(cols[1])
	f.close()

	if DRAW == 'all':
		for chrom in chrDic:
			dic[chrom] = [[0, chrDic[chrom]+1]]
			taille_total += chrDic[chrom]
	else:
		dico = {}
		regions = DRAW.split('-')
		for region in regions:
			elts = region.split(':')
			if not elts[0] in dico:
				dico[elts[0]] = []

			if len(elts) == 3: # we have coordinate
				dico[elts[0]].append([int(elts[1]), int(elts[2])+1])
			else: # draw the whole chromosome
				try:
					dico[elts[0]].append([0, chrDic[elts[0]]+1])
				except Exception, e:
					print (e)
					print ("No chromosome name \""+elts[0]+"\" found in the file : "+CHR+".")

		# to merge overlapping regions
		chrom = ''
		for n in dico:
			list_2_sort = list(dico[n])
			sorted_liste = list(sorted(list_2_sort, key=operator.itemgetter(0)))
			for k in sorted_liste:
				if n in dic:
					if fin < int(k[0]):
						dic[n].append([debut, fin])
						taille_total = taille_total + (fin-debut) + 1
						debut = int(k[0])
						fin = int(k[1])
					else:
						fin = int(k[1])
				else:
					dic[n] = []
					debut = int(k[0])
					fin = int(k[1])
			dic[n].append([debut, fin])
			taille_total = taille_total + (fin-debut) + 1

	if taille_total < 5000:
		print 'Unit : 1 b'
		unit = '10'
	elif taille_total < 50000:
		print 'Unit : 1 Kb'
		unit = '1000'
	elif taille_total < 500000:
		print 'Unit : 10 Kb', taille_total
		unit = '10000'
	elif taille_total < 5000000:
		print 'Unit : 100 Kb'
		unit = '100000'
	elif taille_total < 50000000:
		print 'Unit : 1 Mb'
		unit = '1000000'
	else:
		print 'Unit : 10 Mb'
		unit = '10000000'

	#creation of chromsomes to draw
	chr_order = '^'
	chr_name = ''
	char = 0
	file = open(CHR)
	for line in file:
		data = line.split()
		if data:
			if data[0] in dic:
				list_2_sort = list(dic[data[0]])
				sorted_liste = list(sorted(list_2_sort, key=operator.itemgetter(1)))
				for n in sorted_liste:
					if char > 25:
						mot_char = chr(int(char/25)+96)+chr((char%25)+97)
					else:
						mot_char = chr(char+97)
					chr_order = chr_order+','+str(mot_char)
					if taille_total < 5000:
						if chr_name == '':
							chr_name = data[0]+'['+mot_char+']:%f-%f' % (int(n[0])/float(10), int(n[1])/float(10))
						else:
							chr_name = chr_name+';'+data[0]+'['+mot_char+']:%f-%f' % (int(n[0])/float(10), int(n[1])/float(10))
					elif taille_total < 50000:
						if chr_name == '':
							chr_name = data[0]+'['+mot_char+']:%f-%f' % (int(n[0])/float(1000), int(n[1])/float(1000))
						else:
							chr_name = chr_name+';'+data[0]+'['+mot_char+']:%f-%f' % (int(n[0])/float(1000), int(n[1])/float(1000))
					elif taille_total < 500000:
						if chr_name == '':
							chr_name = data[0]+'['+mot_char+']:%f-%f' % (int(n[0])/float(10000), int(n[1])/float(10000))
						else:
							chr_name = chr_name+';'+data[0]+'['+mot_char+']:%f-%f' % (int(n[0])/float(10000), int(n[1])/float(10000))
					elif taille_total < 5000000:
						if chr_name == '':
							chr_name = data[0]+'['+mot_char+']:%f-%f' % (int(n[0])/float(100000), int(n[1])/float(100000))
						else:
							chr_name = chr_name+';'+data[0]+'['+mot_char+']:%f-%f' % (int(n[0])/float(100000), int(n[1])/float(100000))
					elif taille_total < 50000000:
						if chr_name == '':
							chr_name = data[0]+'['+mot_char+']:%f-%f' % (int(n[0])/float(1000000), int(n[1])/float(1000000))
						else:
							chr_name = chr_name+';'+data[0]+'['+mot_char+']:%f-%f' % (int(n[0])/float(1000000), int(n[1])/float(1000000))
					else:
						if chr_name == '':
							chr_name = data[0]+'['+mot_char+']:%f-%f' % (int(n[0])/float(1000000), int(n[1])/float(1000000))
						else:
							chr_name = chr_name+';'+data[0]+'['+mot_char+']:%f-%f' % (int(n[0])/float(1000000), int(n[1])/float(1000000))
					char = char + 1
	print chr_name
	print chr_order
	return [chr_name, chr_order, unit, dic]

def create_ideogram(FILE, LABEL):
	outfile = open(FILE,'w')
	outfile.write('<ideogram>\n')
	outfile.write('<spacing>\n')
	outfile.write('default = 0.1u\n')
	outfile.write('</spacing>\n')
	outfile.write('thickness        = 50p\n')
	outfile.write('stroke_thickness = 2\n')
	outfile.write('stroke_color     = black\n')
	outfile.write('fill             = yes\n')
	outfile.write('fill_color       = black\n')
	outfile.write('radius         = 0.65r\n')
	if LABEL == 'y':
		outfile.write('show_label     = yes\n')
	else:
		outfile.write('show_label     = no\n')
	outfile.write('label_font     = condensedbold\n')
	outfile.write('label_radius   = dims(ideogram,radius) + 0.35r\n')
	outfile.write('label_size     = 30\n')
	outfile.write('label_parallel = no\n')
	outfile.write('band_stroke_thickness = 0\n')
	outfile.write('show_bands            = yes\n')
	outfile.write('fill_bands            = yes\n')
	outfile.write('band_transparency     = 1\n')
	outfile.write('</ideogram>\n')
	outfile.close()

#create ticks.conf
def create_ticks(FILE, UNIT):
	outfile = open(FILE,'w')
	outfile.write('show_ticks          = yes\n')
	outfile.write('show_tick_labels    = yes\n')
	outfile.write('chrticklabels       = yes\n')
	outfile.write('chrticklabelfont    = default\n')
	outfile.write('grid_start         = dims(ideogram,radius_inner)-0.5r\n')
	outfile.write('grid_end           = dims(ideogram,radius_outer)+100\n')
	outfile.write('<ticks>\n')
	outfile.write('skip_first_label     = no\n')
	outfile.write('skip_last_label      = no\n')
	outfile.write('radius               = dims(ideogram,radius_outer)\n')
	outfile.write('tick_separation      = 2p\n')
	outfile.write('min_label_distance_to_edge = 0p\n')
	outfile.write('label_separation = 5p\n')
	outfile.write('label_offset     = 2p\n')
	outfile.write('label_size = 20p\n')
	if UNIT == '10':
		outfile.write('multiplier = 1\n')
	elif UNIT == '1000':
		outfile.write('multiplier = 1e-3\n')
	elif UNIT == '10000':
		outfile.write('multiplier = 1e-4\n')
	elif UNIT == '100000':
		outfile.write('multiplier = 1e-5\n')
	elif UNIT == '1000000':
		outfile.write('multiplier = 1e-6\n')
	elif UNIT == '10000000':
		outfile.write('multiplier = 1e-7\n')
	outfile.write('color = black\n')
	outfile.write('<tick>\n')
	outfile.write('spacing        = 0.5u\n')
	outfile.write('size           = 5p\n')
	outfile.write('thickness      = 2p\n')
	outfile.write('color          = black\n')
	outfile.write('show_label     = no\n')
	outfile.write('label_size     = 8p\n')
	outfile.write('label_offset   = 0p\n')
	outfile.write('format         = %.0f\n')
	outfile.write('grid           = yes\n')
	outfile.write('grid_color     = grey\n')
	outfile.write('grid_thickness = 1p\n')
	outfile.write('</tick>\n')
	outfile.write('<tick>\n')
	outfile.write('spacing        = 1u\n')
	outfile.write('size           = 8p\n')
	outfile.write('thickness      = 2p\n')
	outfile.write('color          = black\n')
	outfile.write('show_label     = yes\n')
	outfile.write('label_size     = 12p\n')
	outfile.write('label_offset   = 0p\n')
	outfile.write('format         = %.0f\n')
	outfile.write('grid           = yes\n')
	outfile.write('grid_color     = dgrey\n')
	outfile.write('grid_thickness = 1p\n')
	outfile.write('</tick>\n')
	outfile.write('</ticks>\n')
	outfile.close()

def create_sous_tile(TILE, OUT, DIC):
	outfile = open(OUT,'w')
	file = open(TILE)
	for line in file:
		data = line.split()
		if data != []:
			if data[0] in DIC:
				for j in DIC[data[0]]:
					if int(data[1]) >= int(j[0]) and int(data[1]) <= int(j[1]): #le debut est dedans
						if int(data[2]) <= int(j[1]): #la fin est dedans
							outfile.write(line)
						else:# la fin est dehors
							outfile.write(data[0]+'\t'+data[1]+'\t'+str(j[1])+'\t'+data[3]+'\n')
					elif int(data[2]) >= int(j[0]) and int(data[2]) <= int(j[1]): #la fin est dedans mais pas le debut
						outfile.write(data[0]+'\t'+str(j[0])+'\t'+data[2]+'\t'+data[3]+'\n')
					elif int(data[1]) <= int(j[0]) and int(data[2]) >= int(j[1]):
						outfile.write(data[0]+'\t'+str(j[0])+'\t'+str(j[1])+'\t'+data[3]+'\n')
	outfile.close()

def create_sous_object(OBJECT, OUT, DIC):
	outfile = open(OUT,'w')
	file = open(OBJECT)
	for line in file:
		data = line.split()
		if data != []:
			if data[0] in DIC:
				for j in DIC[data[0]]:
					if int(data[2]) >= int(j[0]) and int(data[2]) <= int(j[1]):#la fin est dedans
						outfile.write(line)
					elif int(data[1]) >= int(j[0]) and int(data[1]) <= int(j[1]):#le debut est dedans
						outfile.write(line)
	outfile.close()
	file.close()

def create_kar(KAR, OUT, DIC):
	outfile = open(OUT, 'w')
	file = open(KAR)
	for line in file:
		data = line.split()
		if data != []:
			if data[0] == 'chr':
				if data[2] in DIC:
					outfile.write(line)
			elif data[0] == 'band':
				if data[1] in DIC:
					for j in DIC[data[1]]:
						if int(data[5]) >= int(j[0]) and int(data[5]) <= int(j[1]):#la fin est dedans
							outfile.write(line)
						elif int(data[4]) >= int(j[0]) and int(data[4]) <= int(j[1]):#le debut est dedans
							outfile.write(line)
	outfile.close()
	file.close()

def create_conf(DIC, CONF, IDEOGRAM, TICKS, CONFIG, OUT, CHR, ORDER, UNIT, COV, SCAFF, DISCORD, FRF, FF, RR, INS, DELET, CHR_RR, CHR_FR, CHR_RF, CHR_FF, READ_FR, READ_RF, READ_FF, READ_RR, READ_INS, READ_DELET, READ_CHR_RR, READ_CHR_RF, READ_CHR_FR, READ_CHR_FF, TEXT, KAR):

	liste = []

	configuration = ConfigParser.RawConfigParser()
	configuration.read(CONFIG)
	configuration.get('General','out_kar')

	create_kar(os.path.realpath(configuration.get('General','out_kar')), KAR, DIC)

	outfile = open(CONF,'w')
	outfile.write('<colors>\n')
	outfile.write('<<include etc/colors.conf>>\n')
	outfile.write('</colors>\n')
	outfile.write('<fonts>\n')
	outfile.write('<<include etc/fonts.conf>>\n')
	outfile.write('</fonts>\n')
	outfile.write('<<include %s>>\n' % os.path.realpath(IDEOGRAM))
	outfile.write('<<include %s>>\n' % os.path.realpath(TICKS))
	outfile.write('karyotype = %s\n' % os.path.realpath(KAR))
	outfile.write('<image>\n')
	outfile.write('dir = %s\n' % '/'.join(os.path.realpath(OUT).split('/')[0:-1]))
	outfile.write('file = %s\n' % os.path.realpath(OUT).split('/')[-1])
	outfile.write('24bit = yes\n')
	outfile.write('radius = 1500p\n')
	outfile.write('background = white\n')
	outfile.write('angle_offset = -90\n')
	outfile.write('auto_alpha_colors = yes\n')
	outfile.write('auto_alpha_steps = 5\n')
	outfile.write('</image>\n')
	outfile.write('chromosomes_units = %s\n' % UNIT)
	outfile.write('chromosomes = %s\n' % CHR)
	outfile.write('chromosomes_order = %s\n' % ORDER)

	if READ_FR == 'y' or READ_RF == 'y' or READ_FF == 'y' or READ_RR == 'y' or READ_INS == 'y' or READ_DELET == 'y' or READ_CHR_RR == 'y' or READ_CHR_RF == 'y' or READ_CHR_FR == 'y' or READ_CHR_FF == 'y' or FRF == 'y' or FF == 'y' or RR == 'y' or DELET == 'y' or INS == 'y' or CHR_FF == 'y' or CHR_RF == 'y' or CHR_FR == 'y' or CHR_RR == 'y':
		if (configuration.get('General','read_chr_rr') == 'yes' and READ_CHR_RR == 'y') or (configuration.get('General','read_chr_fr') == 'yes' and READ_CHR_FR == 'y') or (configuration.get('General','read_chr_ff') == 'yes' and READ_CHR_FF == 'y') or (configuration.get('General','read_chr_rf') == 'yes' and READ_CHR_RF == 'y') or (configuration.get('General','read_ins') == 'yes' and READ_INS == 'y') or (configuration.get('General','read_del') == 'yes' and READ_DELET == 'y') or (configuration.get('General','read_rr') == 'yes' and READ_RR == 'y') or (configuration.get('General','read_rf') == 'yes' and READ_RF == 'y') or (configuration.get('General','read_fr') == 'yes' and READ_FR == 'y') or (configuration.get('General','read_ff') == 'yes' and READ_FF == 'y') or (configuration.get('General','chr_rr') == 'yes' and CHR_RR == 'y') or (configuration.get('General','chr_fr') == 'yes' and CHR_FR == 'y') or (configuration.get('General','chr_ff') == 'yes' and CHR_FF == 'y') or (configuration.get('General','chr_rf') == 'yes' and CHR_RF == 'y') or (configuration.get('General','ins') == 'yes' and INS == 'y') or (configuration.get('General','delet') == 'yes' and DELET == 'y') or (configuration.get('General','rr') == 'yes' and RR == 'y') or (configuration.get('General','frf') == 'yes' and FRF == 'y') or (configuration.get('General','ff') == 'yes' and FF == 'y'):
			outfile.write('<links>\n')
			outfile.write('z = 1\n')
			outfile.write('radius = 0.985r\n')
			outfile.write('crest  = 0\n')
			if FRF == 'y':
				if (configuration.get('General','frf')) == 'yes':
					frf = tempfile.NamedTemporaryFile().name
					# print frf
					liste.append(frf)
					cree_link_final(DIC, configuration.get('Discord_link','frf'), frf)
					outfile.write('<link zone_frf_link>\n')
					outfile.write('bezier_radius        = 0.35r\n')
					outfile.write('bezier_radius_purity = 0.5\n')
					outfile.write('show = yes\n')
					outfile.write('color = lblue\n')
					outfile.write('thickness = 1\n')
					outfile.write('file = %s\n' % frf)
					outfile.write('ribbon = yes\n')
					outfile.write('</link>\n')
				else:
					print 'Warning, frf layer has been requested but files are not found in the config file'
			if FF == 'y':
				if (configuration.get('General','ff')) == 'yes':
					ff = tempfile.NamedTemporaryFile().name
					# print ff
					liste.append(ff)
					cree_link_final(DIC, configuration.get('Discord_link','ff'), ff)
					outfile.write('<link zone_ff_link>\n')
					outfile.write('bezier_radius        = 0.5r\n')
					outfile.write('bezier_radius_purity = 0.5\n')
					outfile.write('show = yes\n')
					outfile.write('color = lgreen\n')
					outfile.write('thickness = 1\n')
					outfile.write('file = %s\n' % ff)
					outfile.write('ribbon = yes\n')
					outfile.write('</link>\n')
				else:
					print 'Warning, ff layer has been requested but files are not found in the config file'
			if RR == 'y':
				if (configuration.get('General','rr')) == 'yes':
					rr = tempfile.NamedTemporaryFile().name
					# print rr
					liste.append(rr)
					cree_link_final(DIC, configuration.get('Discord_link','rr'), rr)
					outfile.write('<link zone_rr_link>\n')
					outfile.write('bezier_radius        = 0.45r\n')
					outfile.write('bezier_radius_purity = 0.5\n')
					outfile.write('show = yes\n')
					outfile.write('color = lpurple\n')
					outfile.write('thickness = 1\n')
					outfile.write('file = %s\n' % rr)
					outfile.write('ribbon = yes\n')
					outfile.write('</link>\n')
				else:
					print 'Warning, rr layer has been requested but files are not found in the config file'
			if INS == 'y':
				if (configuration.get('General','ins')) == 'yes':
					ins = tempfile.NamedTemporaryFile().name
					# print ins
					liste.append(ins)
					cree_link_final(DIC, configuration.get('Discord_link','ins'), ins)
					outfile.write('<link zone_ins_link>\n')
					outfile.write('bezier_radius        = 0.3r\n')
					outfile.write('bezier_radius_purity = 0.5\n')
					outfile.write('show = yes\n')
					outfile.write('color = lorange\n')
					outfile.write('thickness = 1\n')
					outfile.write('file = %s\n' % ins)
					outfile.write('ribbon = yes\n')
					outfile.write('</link>\n')
				else:
					print 'Warning, ins layer has been requested but files are not found in the config file'
			if DELET == 'y':
				if (configuration.get('General','delet')) == 'yes':
					delet = tempfile.NamedTemporaryFile().name
					# print delet
					liste.append(delet)
					cree_link_final(DIC, configuration.get('Discord_link','delet'), delet)
					outfile.write('<link zone_delet_link>\n')
					outfile.write('bezier_radius        = 0.4r\n')
					outfile.write('bezier_radius_purity = 0.5\n')
					outfile.write('show = yes\n')
					outfile.write('color = lred\n')
					outfile.write('thickness = 1\n')
					outfile.write('file = %s\n' % delet)
					outfile.write('ribbon = yes\n')
					outfile.write('</link>\n')
				else:
					print 'Warning, delet layer has been requested but files are not found in the config file'
			if CHR_RR == 'y':
				if (configuration.get('General','chr_rr')) == 'yes':
					chr_rr = tempfile.NamedTemporaryFile().name
					# print chr_rr
					liste.append(chr_rr)
					cree_link_final(DIC, configuration.get('Discord_link','chr_rr'), chr_rr)
					outfile.write('<link zone_chr_rr_link>\n')
					outfile.write('bezier_radius        = 0.45r\n')
					outfile.write('bezier_radius_purity = 0.5\n')
					outfile.write('show = yes\n')
					outfile.write('color = lpurple\n')
					outfile.write('thickness = 1\n')
					outfile.write('file = %s\n' % chr_rr)
					outfile.write('ribbon = yes\n')
					outfile.write('</link>\n')
				else:
					print 'Warning, chr_rr layer has been requested but files are not found in the config file'
			if CHR_RF == 'y':
				if (configuration.get('General','chr_rf')) == 'yes':
					chr_rf = tempfile.NamedTemporaryFile().name
					# print chr_rf
					liste.append(chr_rf)
					cree_link_final(DIC, configuration.get('Discord_link','chr_rf'), chr_rf)
					outfile.write('<link zone_chr_rf_link>\n')
					outfile.write('bezier_radius        = 0.4r\n')
					outfile.write('bezier_radius_purity = 0.5\n')
					outfile.write('show = yes\n')
					outfile.write('color = lred\n')
					# if (configuration.get('General','orient')) == 'rf':
						# outfile.write('color = red\n')
					# elif (configuration.get('General','orient')) == 'fr':
						# outfile.write('color = blue\n')
					# else:
						# sys.exit('Unsupported orientation: only fr or rf are supported')
					outfile.write('thickness = 1\n')
					outfile.write('file = %s\n' % chr_rf)
					outfile.write('ribbon = yes\n')
					outfile.write('</link>\n')
				else:
					print 'Warning, chr_rf layer has been requested but files are not found in the config file'
			if CHR_FR == 'y':
				if (configuration.get('General','chr_fr')) == 'yes':
					chr_fr = tempfile.NamedTemporaryFile().name
					# print chr_fr
					liste.append(chr_fr)
					cree_link_final(DIC, configuration.get('Discord_link','chr_fr'), chr_fr)
					outfile.write('<link zone_chr_fr_link>\n')
					outfile.write('bezier_radius        = 0.35r\n')
					outfile.write('bezier_radius_purity = 0.5\n')
					outfile.write('show = yes\n')
					outfile.write('color = lblue\n')
					# if (configuration.get('General','orient')) == 'rf':
						# outfile.write('color = blue\n')
					# elif (configuration.get('General','orient')) == 'fr':
						# outfile.write('color = red\n')
					# else:
						# sys.exit('Unsupported orientation: only fr or rf are supported')
					outfile.write('thickness = 1\n')
					outfile.write('file = %s\n' % chr_fr)
					outfile.write('ribbon = yes\n')
					outfile.write('</link>\n')
				else:
					print 'Warning, chr_fr layer has been requested but files are not found in the config file'
			if CHR_FF == 'y':
				if (configuration.get('General','chr_ff')) == 'yes':
					chr_ff = tempfile.NamedTemporaryFile().name
					# print chr_ff
					liste.append(chr_ff)
					cree_link_final(DIC, configuration.get('Discord_link','chr_ff'), chr_ff)
					outfile.write('<link zone_chr_ff_link>\n')
					outfile.write('bezier_radius        = 0.5r\n')
					outfile.write('bezier_radius_purity = 0.5\n')
					outfile.write('show = yes\n')
					outfile.write('color = lgreen\n')
					outfile.write('thickness = 1\n')
					outfile.write('file = %s\n' % chr_ff)
					outfile.write('ribbon = yes\n')
					outfile.write('</link>\n')
				else:
					print 'Warning, chr_ff layer has been requested but files are not found in the config file'
			if READ_FR == 'y':
				if (configuration.get('General','read_fr')) == 'yes':
					read_fr = tempfile.NamedTemporaryFile().name
					# print read_fr
					liste.append(read_fr)
					cree_link_final(DIC, configuration.get('Read_link','fr'), read_fr)
					outfile.write('<link read_fr_link>\n')
					if (configuration.get('General','orient')) == 'rf':
						outfile.write('perturb        = no\n')
						outfile.write('bezier_radius        = 0.35r\n')
						outfile.write('bezier_radius_purity = 0.5\n')
						outfile.write('show = yes\n')
						outfile.write('color = dblue\n')
						outfile.write('z = 5\n')
					else:
						outfile.write('perturb        = no\n')
						outfile.write('bezier_radius        = 0.55r\n')
						outfile.write('bezier_radius_purity = 0.5\n')
						outfile.write('show = yes\n')
						outfile.write('color = grey\n')
						outfile.write('z = 1\n')
					outfile.write('thickness = 1\n')
					outfile.write('file = %s\n' % read_fr)
					outfile.write('ribbon = yes\n')
					outfile.write('</link>\n')
				else:
					print 'Warning, read_fr layer has been requested but files are not found in the config file'
			if READ_RF == 'y':
				if (configuration.get('General','read_rf')) == 'yes':
					read_rf = tempfile.NamedTemporaryFile().name
					# print read_rf
					liste.append(read_rf)
					cree_link_final(DIC, configuration.get('Read_link','rf'), read_rf)
					outfile.write('<link read_rf_link>\n')
					if (configuration.get('General','orient')) == 'rf':
						outfile.write('perturb        = no\n')
						outfile.write('bezier_radius        = 0.55r\n')
						outfile.write('bezier_radius_purity = 0.5\n')
						outfile.write('show = yes\n')
						outfile.write('color = grey\n')
						outfile.write('z = 1\n')
					else:
						outfile.write('perturb        = no\n')
						outfile.write('bezier_radius        = 0.35r\n')
						outfile.write('bezier_radius_purity = 0.5\n')
						outfile.write('show = yes\n')
						outfile.write('color = dblue\n')
						outfile.write('z = 5\n')
					outfile.write('thickness = 1\n')
					outfile.write('file = %s\n' % read_rf)
					outfile.write('ribbon = yes\n')
					outfile.write('</link>\n')
				else:
					print 'Warning, read_rf layer has been requested but files are not found in the config file'
			if READ_FF == 'y':
				if (configuration.get('General','read_ff')) == 'yes':
					read_ff = tempfile.NamedTemporaryFile().name
					# print read_ff
					liste.append(read_ff)
					cree_link_final(DIC, configuration.get('Read_link','ff'), read_ff)
					outfile.write('<link read_ff_link>\n')
					outfile.write('perturb        = no\n')
					outfile.write('bezier_radius        = 0.5r\n')
					outfile.write('bezier_radius_purity = 0.5\n')
					outfile.write('show = yes\n')
					outfile.write('color = dgreen\n')
					outfile.write('z = 10\n')
					outfile.write('thickness = 1\n')
					outfile.write('file = %s\n' % read_ff)
					outfile.write('ribbon = yes\n')
					outfile.write('</link>\n')
				else:
					print 'Warning, read_ff layer has been requested but files are not found in the config file'
			if READ_RR == 'y':
				if (configuration.get('General','read_rr')) == 'yes':
					read_rr = tempfile.NamedTemporaryFile().name
					# print read_rr
					liste.append(read_rr)
					cree_link_final(DIC, configuration.get('Read_link','rr'), read_rr)
					outfile.write('<link read_rr_link>\n')
					outfile.write('perturb        = no\n')
					outfile.write('bezier_radius        = 0.45r\n')
					outfile.write('bezier_radius_purity = 0.5\n')
					outfile.write('show = yes\n')
					outfile.write('color = dpurple\n')
					outfile.write('z = 10\n')
					outfile.write('thickness = 1\n')
					outfile.write('file = %s\n' % read_rr)
					outfile.write('ribbon = yes\n')
					outfile.write('</link>\n')
				else:
					print 'Warning, read_rr layer has been requested but files are not found in the config file'
			if READ_INS == 'y':
				if (configuration.get('General','read_ins')) == 'yes':
					read_ins = tempfile.NamedTemporaryFile().name
					# print read_ins
					liste.append(read_ins)
					cree_link_final(DIC, configuration.get('Read_link','ins'), read_ins)
					outfile.write('<link read_ins_link>\n')
					outfile.write('perturb        = no\n')
					outfile.write('bezier_radius        = 0.6r\n')
					outfile.write('bezier_radius_purity = 0.5\n')
					outfile.write('show = yes\n')
					outfile.write('color = dorange\n')
					outfile.write('z = 10\n')
					outfile.write('thickness = 1\n')
					outfile.write('file = %s\n' % read_ins)
					outfile.write('ribbon = yes\n')
					outfile.write('</link>\n')
				else:
					print 'Warning, read_ins layer has been requested but files are not found in the config file'
			if READ_DELET == 'y':
				if (configuration.get('General','read_del')) == 'yes':
					read_del = tempfile.NamedTemporaryFile().name
					# print read_del
					liste.append(read_del)
					cree_link_final(DIC, configuration.get('Read_link','del'), read_del)
					outfile.write('<link read_delet_link>\n')
					outfile.write('perturb        = no\n')
					outfile.write('bezier_radius        = 0.4r\n')
					outfile.write('bezier_radius_purity = 0.5\n')
					outfile.write('show = yes\n')
					outfile.write('color = dred\n')
					outfile.write('z = 10\n')
					outfile.write('thickness = 1\n')
					outfile.write('file = %s\n' % read_del)
					outfile.write('ribbon = yes\n')
					outfile.write('</link>\n')
				else:
					print 'Warning, read_del layer has been requested but files are not found in the config file'
			if READ_CHR_FR == 'y':
				if (configuration.get('General','read_chr_fr')) == 'yes':
					read_chr_fr = tempfile.NamedTemporaryFile().name
					# print read_chr_fr
					liste.append(read_chr_fr)
					cree_link_final(DIC, configuration.get('Read_link','chr_fr'), read_chr_fr)
					outfile.write('<link read_chr_fr_link>\n')
					outfile.write('perturb        = no\n')
					outfile.write('bezier_radius        = 0.35r\n')
					outfile.write('bezier_radius_purity = 0.5\n')
					outfile.write('show = yes\n')
					outfile.write('color = dblue\n')
					# if (configuration.get('General','orient')) == 'rf':
						# outfile.write('color = dblue\n')
					# elif (configuration.get('General','orient')) == 'fr':
						# outfile.write('color = dred\n')
					# else:
						# sys.exit('Unsupported orientation: only fr or rf are supported')
					outfile.write('z = 5\n')
					outfile.write('thickness = 1\n')
					outfile.write('file = %s\n' % read_chr_fr)
					outfile.write('ribbon = yes\n')
					outfile.write('</link>\n')
				else:
					print 'Warning, read_chr_fr layer has been requested but files are not found in the config file'
			if READ_CHR_RF == 'y':
				if (configuration.get('General','read_chr_rf')) == 'yes':
					read_chr_rf = tempfile.NamedTemporaryFile().name
					# print read_chr_rf
					liste.append(read_chr_rf)
					cree_link_final(DIC, configuration.get('Read_link','chr_rf'), read_chr_rf)
					outfile.write('<link read_chr_rf_link>\n')
					outfile.write('perturb        = no\n')
					outfile.write('bezier_radius        = 0.4r\n')
					outfile.write('bezier_radius_purity = 0.5\n')
					outfile.write('show = yes\n')
					outfile.write('color = dred\n')
					# if (configuration.get('General','orient')) == 'rf':
						# outfile.write('color = dred\n')
					# elif (configuration.get('General','orient')) == 'fr':
						# outfile.write('color = dblue\n')
					# else:
						# sys.exit('Unsupported orientation: only fr or rf are supported')
					outfile.write('z = 5\n')
					outfile.write('thickness = 1\n')
					outfile.write('file = %s\n' % read_chr_rf)
					outfile.write('ribbon = yes\n')
					outfile.write('</link>\n')
				else:
					print 'Warning, read_chr_rf layer has been requested but files are not found in the config file'
			if READ_CHR_FF == 'y':
				if (configuration.get('General','read_chr_ff')) == 'yes':
					read_chr_ff = tempfile.NamedTemporaryFile().name
					# print read_chr_ff
					liste.append(read_chr_ff)
					cree_link_final(DIC, configuration.get('Read_link','chr_ff'), read_chr_ff)
					outfile.write('<link read_chr_ff_link>\n')
					outfile.write('perturb        = no\n')
					outfile.write('bezier_radius        = 0.5r\n')
					outfile.write('bezier_radius_purity = 0.5\n')
					outfile.write('show = yes\n')
					outfile.write('color = dgreen\n')
					outfile.write('z = 10\n')
					outfile.write('thickness = 1\n')
					outfile.write('file = %s\n' % read_chr_ff)
					outfile.write('ribbon = yes\n')
					outfile.write('</link>\n')
				else:
					print 'Warning, read_chr_ff layer has been requested but files are not found in the config file'
			if READ_CHR_RR == 'y':
				if (configuration.get('General','read_chr_rr')) == 'yes':
					read_chr_rr = tempfile.NamedTemporaryFile().name
					# print read_chr_rr
					liste.append(read_chr_rr)
					cree_link_final(DIC, configuration.get('Read_link','chr_rr'), read_chr_rr)
					outfile.write('<link read_chr_rr_link>\n')
					outfile.write('perturb        = no\n')
					outfile.write('bezier_radius        = 0.45r\n')
					outfile.write('bezier_radius_purity = 0.5\n')
					outfile.write('show = yes\n')
					outfile.write('color = dpurple\n')
					outfile.write('z = 10\n')
					outfile.write('thickness = 1\n')
					outfile.write('file = %s\n' % read_chr_rr)
					outfile.write('ribbon = yes\n')
					outfile.write('</link>\n')
				else:
					print 'Warning, read_chr_rr layer has been requested but files are not found in the config file'
			outfile.write('</links>\n')

	if COV == 'y' or SCAFF == 'y' or DISCORD == 'y' or TEXT == 'y':
		if TEXT == 'y' or (configuration.get('General','prop') == 'yes' and DISCORD == 'y') or (configuration.get('General','cov') == 'yes' and COV == 'y') or (configuration.get('General','scaff_tile') == 'yes' and SCAFF == 'y'):
			outfile.write('<plots>\n')
			if TEXT == 'y':
				text = tempfile.NamedTemporaryFile().name
				# print text
				liste.append(text)
				cree_text_final(DIC, os.path.realpath(configuration.get('General','out_n')), text)
				outfile.write('<plot>\n')
				outfile.write('z = 1\n')
				outfile.write('show             = yes\n')
				# outfile.write('label_rotate     = no\n')
				outfile.write('type             = text\n')
				outfile.write('color            = black\n')
				outfile.write('file             = %s\n' % text)
				outfile.write('r0               = 1r\n')
				outfile.write('r1               = 2r\n')
				outfile.write('label_size       = 20p\n')
				outfile.write('padding          = 10p\n')
				outfile.write('rpadding         = 0p\n')
				outfile.write('show_links     = yes\n')
				outfile.write('link_dims      = 0p,180p,20p,20p,10p\n')
				outfile.write('link_thickness = 2p\n')
				outfile.write('link_color     = red\n')
				outfile.write('label_snuggle         = yes\n')
				outfile.write('max_snuggle_distance  = 3r\n')
				outfile.write('snuggle_tolerance     = 0.25r\n')
				outfile.write('snuggle_sampling      = 2\n')
				outfile.write('snuggle_link_overlap_test = yes\n')
				outfile.write('snuggle_link_overlap_tolerance = 2p\n')
				outfile.write('snuggle_refine        = yes\n')
				outfile.write('label_rotate = yes\n')
				outfile.write('</plot>\n')
			if DISCORD == 'y':
				if configuration.get('General','prop') == 'yes':
					prop_file = tempfile.NamedTemporaryFile().name
					# print prop_file
					liste.append(prop_file)
					create_sous_object(os.path.realpath(configuration.get('Proportion','prop')), prop_file, DIC)
					outfile.write('<plot>\n')
					outfile.write('z = 1\n')
					outfile.write('show  = yes\n')
					outfile.write('type  = scatter\n')
					outfile.write('file  = %s\n' % prop_file)
					outfile.write('r1    = 1.1r\n')
					outfile.write('r0    = 1.04r\n')
					outfile.write('max   = 1\n')
					outfile.write('min   = -0.2\n')
					outfile.write('axis           = yes\n')
					outfile.write('axis_color     = lgrey\n')
					outfile.write('axis_thickness = 1\n')
					outfile.write('axis_spacing   = 0.25\n')
					outfile.write('glyph            = circle\n')
					outfile.write('glyph_size       = 5\n')
					outfile.write('color            = red\n')
					outfile.write('fill_color       = red\n')
					outfile.write('stroke_color     = red\n')
					outfile.write('stroke_thickness = 0\n')
					outfile.write('<rules>\n')
					outfile.write('<rule>\n')
					outfile.write('importance = 100\n')
					outfile.write('condition    = _VALUE_ < 0\n')
					outfile.write('color        = black\n')
					outfile.write('fill_color   = black\n')
					outfile.write('glyph        = circle\n')
					outfile.write('glyph_size   = 5\n')
					outfile.write('</rule>\n')
					outfile.write('<rule>\n')
					outfile.write('importance = 90\n')
					outfile.write('condition    = _VALUE_ < 0.25\n')
					outfile.write('color        = green\n')
					outfile.write('fill_color   = green\n')
					outfile.write('glyph        = circle\n')
					outfile.write('glyph_size   = 5\n')
					outfile.write('</rule>\n')
					outfile.write('<rule>\n')
					outfile.write('importance = 80\n')
					outfile.write('condition    = _VALUE_ < 0.50\n')
					outfile.write('color        = yellow\n')
					outfile.write('fill_color   = yellow\n')
					outfile.write('glyph        = circle\n')
					outfile.write('glyph_size   = 5\n')
					outfile.write('</rule>\n')
					outfile.write('<rule>\n')
					outfile.write('importance = 80\n')
					outfile.write('condition    = _VALUE_ < 0.75\n')
					outfile.write('color        = orange\n')
					outfile.write('fill_color   = orange\n')
					outfile.write('glyph        = circle\n')
					outfile.write('glyph_size   = 5\n')
					outfile.write('</rule>\n')
					outfile.write('</rules>\n')
					outfile.write('</plot>\n')
				else:
					print 'Warning, discordant read proportion layer has been requested but files are not found in the config file'
			if COV == 'y':
				if configuration.get('General','cov') == 'yes':
					cov_file = tempfile.NamedTemporaryFile().name
					# print cov_file
					liste.append(cov_file)
					create_sous_object(os.path.realpath(configuration.get('Coverage','cov')), cov_file, DIC)
					outfile.write('<plot>\n')
					outfile.write('z = 1\n')
					outfile.write('show  = yes\n')
					outfile.write('type  = scatter\n')
					outfile.write('file  = %s\n' % cov_file)
					outfile.write('r0    = 1.14r\n')
					outfile.write('r1    = 1.28r\n')
					outfile.write('max   = %s\n' % str(configuration.getfloat('Coverage','mean_cov')+(3*configuration.getfloat('Coverage','median_cov'))))
					outfile.write('min   = %s\n' % str(configuration.getfloat('Coverage','mean_cov')-configuration.getfloat('Coverage','median_cov')))
					outfile.write('glyph            = circle\n')
					outfile.write('glyph_size       = 5\n')
					outfile.write('color            = green\n')
					outfile.write('fill_color       = green\n')
					outfile.write('stroke_color     = black\n')
					outfile.write('stroke_thickness = 0\n')
					outfile.write('axis           = yes\n')
					outfile.write('axis_color     = black\n')
					outfile.write('axis_thickness = 1\n')
					outfile.write('axis_spacing   = %s\n' % str(configuration.getfloat('Coverage','median_cov')*4))
					outfile.write('<rules>\n')
					outfile.write('<rule>\n')
					outfile.write('importance = 100\n')
					outfile.write('condition    = _VALUE_ > %s\n' % str(configuration.getfloat('Coverage','mean_cov')+configuration.getfloat('Coverage','median_cov')))
					outfile.write('color        = red\n')
					outfile.write('fill_color   = red\n')
					outfile.write('glyph        = circle\n')
					outfile.write('glyph_size   = 5\n')
					outfile.write('</rule>\n')
					outfile.write('<rule>\n')
					outfile.write('importance = 90\n')
					outfile.write('condition    = _VALUE_ > %s\n' % str(configuration.getfloat('Coverage','mean_cov')+(configuration.getfloat('Coverage','median_cov')/2)))
					outfile.write('color        = orange\n')
					outfile.write('fill_color   = orange\n')
					outfile.write('glyph        = circle\n')
					outfile.write('glyph_size   = 5\n')
					outfile.write('</rule>\n')
					outfile.write('<rule>\n')
					outfile.write('importance = 80\n')
					outfile.write('condition    = _VALUE_ > %s\n' % str(configuration.getfloat('Coverage','mean_cov')+(configuration.getfloat('Coverage','median_cov')/3)))
					outfile.write('color        = yellow\n')
					outfile.write('fill_color   = yellow\n')
					outfile.write('glyph        = circle\n')
					outfile.write('glyph_size   = 5\n')
					outfile.write('</rule>\n')
					outfile.write('<rule>\n')
					outfile.write('importance = 70\n')
					outfile.write('condition    = _VALUE_ > %s\n' % str(configuration.getfloat('Coverage','mean_cov')+(configuration.getfloat('Coverage','median_cov')/4)))
					outfile.write('color        = lgreen\n')
					outfile.write('fill_color   = lgreen\n')
					outfile.write('glyph        = circle\n')
					outfile.write('glyph_size   = 5\n')
					outfile.write('</rule>\n')
					outfile.write('<rule>\n')
					outfile.write('importance = 100\n')
					outfile.write('condition    = _VALUE_ < %s\n' % str(configuration.getfloat('Coverage','mean_cov')-configuration.getfloat('Coverage','median_cov')+0.001))
					outfile.write('color        = black\n')
					outfile.write('fill_color   = black\n')
					outfile.write('glyph        = circle\n')
					outfile.write('glyph_size   = 5\n')
					outfile.write('</rule>\n')
					outfile.write('<rule>\n')
					outfile.write('importance = 90\n')
					outfile.write('condition    = _VALUE_ < %s\n' % str(configuration.getfloat('Coverage','mean_cov')-(configuration.getfloat('Coverage','median_cov')/2)))
					outfile.write('color        = dblue\n')
					outfile.write('fill_color   = dblue\n')
					outfile.write('glyph        = circle\n')
					outfile.write('glyph_size   = 5\n')
					outfile.write('</rule>\n')
					outfile.write('<rule>\n')
					outfile.write('importance = 80\n')
					outfile.write('condition    = _VALUE_ < %s\n' % str(configuration.getfloat('Coverage','mean_cov')-(configuration.getfloat('Coverage','median_cov')/3)))
					outfile.write('color        = blue\n')
					outfile.write('fill_color   = blue\n')
					outfile.write('glyph        = circle\n')
					outfile.write('glyph_size   = 5\n')
					outfile.write('</rule>\n')
					outfile.write('<rule>\n')
					outfile.write('importance = 70\n')
					outfile.write('condition    = _VALUE_ < %s\n' % str(configuration.getfloat('Coverage','mean_cov')-(configuration.getfloat('Coverage','median_cov')/4)))
					outfile.write('color        = dgreen\n')
					outfile.write('fill_color   = dgreen\n')
					outfile.write('glyph        = circle\n')
					outfile.write('glyph_size   = 5\n')
					outfile.write('</rule>\n')
					outfile.write('</rules>\n')
					outfile.write('</plot>\n')
					#for coverage belong 3*median
					outfile.write('<plot>\n')
					outfile.write('type = histogram\n')
					outfile.write('file  = %s\n' % cov_file)
					outfile.write('orientation = out\n')
					outfile.write('z = 0\n')
					outfile.write('r0    = 1.11r\n')
					outfile.write('r1    = 1.25r\n')
					outfile.write('max   = %s\n' % str(configuration.getfloat('Coverage','mean_cov')+(3*configuration.getfloat('Coverage','median_cov'))))
					outfile.write('min   = %s\n' % str(configuration.getfloat('Coverage','mean_cov')-configuration.getfloat('Coverage','median_cov')))
					outfile.write('color = white\n')
					outfile.write('fill_under = no\n')
					outfile.write('thickness = 1\n')
					outfile.write('extend_bin = no\n')
					outfile.write('axis = no\n')
					outfile.write('background = no\n')
					outfile.write('<rules>\n')
					outfile.write('<rule>\n')
					outfile.write('importance = 100\n')
					outfile.write('condition  = _VALUE_ > %s\n' % str(configuration.getfloat('Coverage','mean_cov')+(3*configuration.getfloat('Coverage','median_cov'))))
					outfile.write('fill_under = yes\n')
					outfile.write('fill_color = red\n')
					outfile.write('color = red\n')
					outfile.write('</rule>\n')
					outfile.write('</rules>\n')
					outfile.write('</plot>\n')
				else:
					print 'Warning, coverage layer has been requested but files are not found in the config file'
			if SCAFF == 'y':
				if configuration.get('General','scaff_tile') == 'yes':
					tile = tempfile.NamedTemporaryFile().name
					# print tile
					liste.append(tile)
					create_sous_tile(os.path.realpath(configuration.get('Scaffold','scaff_tile')), tile, DIC)
					outfile.write('<plot>\n')
					outfile.write('file    = %s\n' % tile)
					outfile.write('show    = yes\n')
					outfile.write('type    = tile\n')
					outfile.write('layers = 15\n')
					outfile.write('margin = 0u\n')
					outfile.write('thickness = 10\n')
					outfile.write('padding = 8\n')
					outfile.write('orientation = out\n')
					outfile.write('stroke_thickness = 1\n')
					outfile.write('stroke_color     = black\n')
					outfile.write('r0    = 0.987r\n')
					outfile.write('r1    = 1r\n')
					outfile.write('background       = no\n')
					outfile.write('</plot>\n')
				else:
					print 'Warning, scaffold layer has been requested but files are not found in the config file'
			outfile.write('</plots>\n')

	outfile.write('anglestep = 0.5\n')
	outfile.write('minslicestep = 10\n')
	outfile.write('beziersamples = 40\n')
	outfile.write('debug = no\n')
	outfile.write('warnings = no\n')
	outfile.write('imagemap = no\n')
	outfile.write('units_ok = bupr\n')
	outfile.write('units_nounit = n\n')
	outfile.close()
	return liste

def cree_text_final(DIC, TEXT, OUT):
	file = open(TEXT)
	outfile = open(OUT,'w')
	for line in file:
		data = line.split()
		if data != []:
			if data[0] in DIC:
				dedans = 0
				for j in DIC[data[0]]:
					if int(data[1]) >= int(j[0]) and int(data[1]) <= int(j[1]):
						outfile.write(line)
						break
	outfile.close()
	file.close()

def cree_link_final(DIC, LINK, OUT):
	file = open(LINK)
	outfile = open(OUT,'w')
	boucle = 1
	while boucle:
		line = file.readline()
		line2 = file.readline()
		data = line.split()
		data2 = line2.split()
		if data != []:
			if data[1] in DIC and data2[1] in DIC:
				dedans = 0
				for j in DIC[data[1]]:
					if int(data[2]) >= int(j[0]) and int(data[3]) <= int(j[1]):
						dedans = 1
						break
				for j in DIC[data2[1]]:
					if int(data2[2]) >= int(j[0]) and int(data2[3]) <= int(j[1]) and dedans == 1:
						dedans = 2
						break
				if dedans == 2:
					outfile.write(line)
					outfile.write(line2)
		else:
			boucle = 0
	outfile.close()
	file.close()


def __main__():
	#Parse Command Line
	parser = optparse.OptionParser(usage="python %prog [options]\n\nProgram designed by Guillaume MARTIN : guillaume.martin@cirad.fr\n\n"
	"This program generate the circos picture")

	# Wrapper options.
	#For input
	parser.add_option( '', '--config', dest='config', default='not_filled', help='The conf file generated by conf4circos.py')
	parser.add_option( '', '--draw', dest='draw', default='all', help='A list of chromosome separated with \'-\'. The position on the chromosome are separated by \':\'. (ex:chr01-chr02:1000000:2000000-chr8), [default: %default]')
	parser.add_option( '', '--cov', dest='cov', default='y', help='Draw coverage plot (y or n), [default: %default]')
	parser.add_option( '', '--scaff', dest='scaff', default='y', help='Draw scaffold position (y or n), [default: %default]')
	parser.add_option( '', '--discord', dest='discord', default='y', help='Draw discordant proportion plot (y or n), [default: %default]')
	parser.add_option( '', '--frf', dest='frf', default='y', help='Draw discordant frf link (y or n), [default: %default]')
	parser.add_option( '', '--ff', dest='ff', default='y', help='Draw discordant ff link (y or n), [default: %default]')
	parser.add_option( '', '--rr', dest='rr', default='y', help='Draw discordant rr link (y or n), [default: %default]')
	parser.add_option( '', '--ins', dest='ins', default='y', help='Draw discordant ins link (y or n), [default: %default]')
	parser.add_option( '', '--delet', dest='delet', default='y', help='Draw discordant delet link (y or n), [default: %default]')
	parser.add_option( '', '--chr_rr', dest='chr_rr', default='y', help='Draw discordant chr_rr link (y or n), [default: %default]')
	parser.add_option( '', '--chr_rf', dest='chr_rf', default='y', help='Draw discordant chr_rf link (y or n), [default: %default]')
	parser.add_option( '', '--chr_fr', dest='chr_fr', default='y', help='Draw discordant chr_fr link (y or n), [default: %default]')
	parser.add_option( '', '--chr_ff', dest='chr_ff', default='y', help='Draw discordant chr_ff link (y or n), [default: %default]')


	parser.add_option( '', '--read_fr', dest='read_fr', default='y', help='Draw read fr link (y or n), [default: %default]')
	parser.add_option( '', '--read_rf', dest='read_rf', default='y', help='Draw read rf link (y or n), [default: %default]')
	parser.add_option( '', '--read_ff', dest='read_ff', default='y', help='Draw read ff link (y or n), [default: %default]')
	parser.add_option( '', '--read_rr', dest='read_rr', default='y', help='Draw read rr link (y or n), [default: %default]')
	parser.add_option( '', '--read_ins', dest='read_ins', default='y', help='Draw read ins link (y or n), [default: %default]')
	parser.add_option( '', '--read_delet', dest='read_delet', default='y', help='Draw read delet link (y or n), [default: %default]')
	parser.add_option( '', '--read_chr_rr', dest='read_chr_rr', default='y', help='Draw read chr_rr link (y or n), [default: %default]')
	parser.add_option( '', '--read_chr_rf', dest='read_chr_rf', default='y', help='Draw read chr_rf link (y or n), [default: %default]')
	parser.add_option( '', '--read_chr_fr', dest='read_chr_fr', default='y', help='Draw read chr_fr link (y or n), [default: %default]')
	parser.add_option( '', '--read_chr_ff', dest='read_chr_ff', default='y', help='Draw read chr_ff link (y or n), [default: %default]')

	parser.add_option( '', '--text', dest='text', default='y', help='Locate N regions, [default: %default]')

	parser.add_option( '', '--labels', dest='labels', default='y', help='Draw reference sequence name, [default: %default]')

	parser.add_option( '', '--out', dest='out', default='circos.png', help='The output file name')
	(options, args) = parser.parse_args()



	pathname = os.path.dirname(sys.argv[0])

	loca_programs = ConfigParser.RawConfigParser()
	loca_programs.read(pathname+'/loca_programs.conf')



	if options.config == 'not_filled':
		sys.exit('--config argument is missing')

	t0 = datetime.datetime.now()
	config = ConfigParser.RawConfigParser()
	config.read(options.config)

	chr_info = define_regions(options.draw, config.get('General','chr'))


	ideogram = tempfile.NamedTemporaryFile().name+'.ideo'
	ticks = tempfile.NamedTemporaryFile().name+'.ticks'
	conf = tempfile.NamedTemporaryFile().name+'.conf'
	kar = tempfile.NamedTemporaryFile().name+'.kar'

	create_ideogram(ideogram, options.labels)
	create_ticks(ticks, chr_info[2])

	liste_rm = create_conf(chr_info[3], conf, ideogram, ticks, options.config, os.path.splitext(options.out)[0], chr_info[0], chr_info[1], chr_info[2], options.cov, options.scaff, options.discord, options.frf, options.ff, options.rr, options.ins, options.delet, options.chr_rr, options.chr_fr, options.chr_rf, options.chr_ff, options.read_fr, options.read_rf, options.read_ff, options.read_rr, options.read_ins, options.read_delet, options.read_chr_rr, options.read_chr_rf, options.read_chr_fr, options.read_chr_ff, options.text, kar)

	run_circos = '%s %s -conf %s' % (loca_programs.get('Programs','perl'), loca_programs.get('Programs','circos'), conf)
	run_job(run_circos, 'Bug in run circos')

	os.system('rm '+ideogram)
	os.system('rm '+ticks)
	os.system('rm '+conf)
	os.system('rm '+kar)
	for n in liste_rm:
		os.system('rm '+n)
	print os.path.splitext(options.out)[0]
	print os.path.splitext(options.out)[0]+'.png'
	print options.out
	if os.path.splitext(options.out)[0]+'.png' != options.out:
		os.system('mv '+os.path.splitext(options.out)[0]+'.png '+options.out)
	print datetime.datetime.now() - t0

if __name__ == "__main__": __main__()
