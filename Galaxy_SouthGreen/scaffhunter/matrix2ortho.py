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

import optparse, os, shutil, subprocess, sys, tempfile, fileinput, ConfigParser, operator, time

def stop_err( msg ):
    sys.stderr.write( "%s\n" % msg )
    sys.exit()

def run_job (cmd_line, ERROR):
	print cmd_line
	try:
		tmp = tempfile.NamedTemporaryFile().name
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


def min_max(FILE):
	file = open(FILE)
	line = file.readline()
	line = file.readline()
	data = line.split()
	val = map(float,data[1:])
	try:
		while True:	val.remove(999999999)
	except: pass
	mini = min(val)
	maxi = max(val)
	i = 1
	for line in file:
		i += 1
		data = line.split()
		if data:
			val = map(float,data[1:i])
			try:
				while True:	val.remove(999999999)
			except: pass
			if val:
				min_val = min(val)
				max_val = max(val)
				if mini > min_val:
					mini = min_val
				if maxi < max_val:
					maxi = max_val
	return [mini, maxi]

def rename(NOM,LISTE):
	if len(str(LISTE.index(NOM))) == 1:
		return '00'+str(LISTE.index(NOM))+'-'+NOM
	elif len(str(LISTE.index(NOM))) == 2:
		return '0'+str(LISTE.index(NOM))+'-'+NOM
	else:
		return str(LISTE.index(NOM))+'-'+NOM

def format_ortho(FILE, DICO, LISTE, LISTE_MARK, MIN, MAX, TYPE, OUT):
	file = open(FILE)
	i = 1
	V1 = (MAX - (1*(MAX/10)))
	V2 = (MAX - (2*(MAX/10)))
	V3 = (MAX - (3*(MAX/10)))
	V4 = (MAX - (4*(MAX/10)))
	V5 = (MAX - (5*(MAX/10)))
	V6 = (MAX - (6*(MAX/10)))
	V7 = (MAX - (7*(MAX/10)))
	V8 = (MAX - (8*(MAX/10)))
	V9 = (MAX - (9*(MAX/10)))
	VV1 = (MIN - (1*(MIN/10)))
	VV2 = (MIN - (2*(MIN/10)))
	VV3 = (MIN - (3*(MIN/10)))
	VV4 = (MIN - (4*(MIN/10)))
	VV5 = (MIN - (5*(MIN/10)))
	VV6 = (MIN - (6*(MIN/10)))
	VV7 = (MIN - (7*(MIN/10)))
	VV8 = (MIN - (8*(MIN/10)))
	VV9 = (MIN - (9*(MIN/10)))
	VVV1 = (9*(MAX/20))
	VVV2 = (8*(MAX/20))
	VVV3 = (7*(MAX/20))
	VVV4 = (6*(MAX/20))
	VVV5 = (5*(MAX/20))
	VVV6 = (4*(MAX/20))
	VVV7 = (3*(MAX/20))
	VVV8 = (2*(MAX/20))
	VVV9 = (1*(MAX/20))
	os.system('echo "legend of the dot-plot"')
	if TYPE == 'REC':
		os.system('echo "[0;'+str(VVV9)+'] black"')
		os.system('echo "]'+str(VVV9)+';'+str(VVV8)+'] red"')
		os.system('echo "]'+str(VVV8)+';'+str(VVV7)+'] pink"')
		os.system('echo "]'+str(VVV7)+';'+str(VVV6)+'] orange"')
		os.system('echo "]'+str(VVV6)+';'+str(VVV5)+'] yellow"')
		os.system('echo "]'+str(VVV5)+';'+str(VVV4)+'] lightgreen"')
		os.system('echo "]'+str(VVV4)+';'+str(VVV3)+'] green"')
		os.system('echo "]'+str(VVV3)+';'+str(VVV2)+'] lightskyblue"')
		os.system('echo "]'+str(VVV2)+';'+str(VVV1)+'] blue3"')
	elif TYPE == 'LOD':
		os.system('echo "['+str(V1)+';'+str(MAX)+'] black"')
		os.system('echo "['+str(V2)+';'+str(V1)+'[ red"')
		os.system('echo "['+str(V3)+';'+str(V2)+'[ pink"')
		os.system('echo "['+str(V4)+';'+str(V3)+'[ orange"')
		os.system('echo "['+str(V5)+';'+str(V4)+'[ yellow"')
		os.system('echo "['+str(V6)+';'+str(V5)+'[ lightgreen"')
		os.system('echo "['+str(V7)+';'+str(V6)+'[ green"')
		os.system('echo "['+str(V8)+';'+str(V7)+'[ lightskyblue"')
		os.system('echo "['+str(V9)+';'+str(V8)+'[ blue3"')
	elif TYPE == 'COR':
		os.system('echo "under the diag"')
		os.system('echo "['+str(V1)+';'+str(MAX)+'] black"')
		os.system('echo "['+str(V2)+';'+str(V1)+'[ red"')
		os.system('echo "['+str(V3)+';'+str(V2)+'[ pink"')
		os.system('echo "['+str(V4)+';'+str(V3)+'[ orange"')
		os.system('echo "['+str(V5)+';'+str(V4)+'[ yellow"')
		os.system('echo "['+str(V6)+';'+str(V5)+'[ lightgreen"')
		os.system('echo "['+str(V7)+';'+str(V6)+'[ green"')
		os.system('echo "['+str(V8)+';'+str(V7)+'[ lightskyblue"')
		os.system('echo "['+str(V9)+';'+str(V8)+'[ blue3"')
		os.system('echo "over the diag"')
		os.system('echo "['+str(MIN)+';'+str(VV1)+'] black"')
		os.system('echo "['+str(VV1)+';'+str(VV2)+'[ red"')
		os.system('echo "['+str(VV2)+';'+str(VV3)+'[ pink"')
		os.system('echo "['+str(VV3)+';'+str(VV4)+'[ orange"')
		os.system('echo "['+str(VV4)+';'+str(VV5)+'[ yellow"')
		os.system('echo "['+str(VV5)+';'+str(VV6)+'[ lightgreen"')
		os.system('echo "['+str(VV6)+';'+str(VV7)+'[ green"')
		os.system('echo "['+str(VV7)+';'+str(VV8)+'[ lightskyblue"')
		os.system('echo "['+str(VV8)+';'+str(VV9)+'[ blue3"')
		
	line = file.readline()
	data = line.split()
	liste_mat_id = data[1:]
	for line in file:
		data = line.split()
		if data:
			if data[0] in LISTE_MARK:
				marker_index = LISTE_MARK.index(data[0])
				nom1 = rename(DICO[data[0]],LISTE)
				pos_mark1d = str((1+marker_index)*10)
				pos_mark1f = str(((1+marker_index)*10)+9)
				# print i
				j = 0
				while j < i:#on ne lit que la diagonale inferieure de la matrice (diagonale comprise)
					j += 1
					name_marker2 = liste_mat_id[j-1]
					if data[j] == '-':
						print 'Warning. Couple  '+data[0]+' '+name_marker2+' has no value'
					elif name_marker2 in LISTE_MARK:
						marker2_index = LISTE_MARK.index(name_marker2)
						valeur = float(data[j])
						if valeur != 999999999:
						# if j == i:
							# print valeur, i, j
							# sys.exit()
							if marker_index <= marker2_index:#meme region ou region ligne avant region colonne
								if TYPE == 'LOD' or TYPE == 'COR':
									if valeur >= V9:
										if valeur >= V1:
											OUT.write('abscisse '+nom1+' '+pos_mark1d+' '+pos_mark1f+' ordonnee '+rename(DICO[name_marker2],LISTE)+' '+str((1+marker2_index)*10)+' '+str(((1+marker2_index)*10)+9)+' black\n')
										elif valeur >= V2:
											OUT.write('abscisse '+nom1+' '+pos_mark1d+' '+pos_mark1f+' ordonnee '+rename(DICO[name_marker2],LISTE)+' '+str((1+marker2_index)*10)+' '+str(((1+marker2_index)*10)+9)+' red\n')
										elif valeur >= V3:
											OUT.write('abscisse '+nom1+' '+pos_mark1d+' '+pos_mark1f+' ordonnee '+rename(DICO[name_marker2],LISTE)+' '+str((1+marker2_index)*10)+' '+str(((1+marker2_index)*10)+9)+' pink\n')
										elif valeur >= V4:
											OUT.write('abscisse '+nom1+' '+pos_mark1d+' '+pos_mark1f+' ordonnee '+rename(DICO[name_marker2],LISTE)+' '+str((1+marker2_index)*10)+' '+str(((1+marker2_index)*10)+9)+' orange\n')
										elif valeur >= V5:
											OUT.write('abscisse '+nom1+' '+pos_mark1d+' '+pos_mark1f+' ordonnee '+rename(DICO[name_marker2],LISTE)+' '+str((1+marker2_index)*10)+' '+str(((1+marker2_index)*10)+9)+' yellow\n')
										elif valeur >= V6:
											OUT.write('abscisse '+nom1+' '+pos_mark1d+' '+pos_mark1f+' ordonnee '+rename(DICO[name_marker2],LISTE)+' '+str((1+marker2_index)*10)+' '+str(((1+marker2_index)*10)+9)+' lightgreen\n')
										elif valeur >= V7:
											OUT.write('abscisse '+nom1+' '+pos_mark1d+' '+pos_mark1f+' ordonnee '+rename(DICO[name_marker2],LISTE)+' '+str((1+marker2_index)*10)+' '+str(((1+marker2_index)*10)+9)+' green\n')
										elif valeur >= V8:
											OUT.write('abscisse '+nom1+' '+pos_mark1d+' '+pos_mark1f+' ordonnee '+rename(DICO[name_marker2],LISTE)+' '+str((1+marker2_index)*10)+' '+str(((1+marker2_index)*10)+9)+' lightskyblue\n')
										else :
											OUT.write('abscisse '+nom1+' '+pos_mark1d+' '+pos_mark1f+' ordonnee '+rename(DICO[name_marker2],LISTE)+' '+str((1+marker2_index)*10)+' '+str(((1+marker2_index)*10)+9)+' blue3\n')
								elif TYPE == 'REC':
									if valeur <= VVV1:
										if valeur <= VVV9:
											OUT.write('abscisse '+nom1+' '+pos_mark1d+' '+pos_mark1f+' ordonnee '+rename(DICO[name_marker2],LISTE)+' '+str((1+marker2_index)*10)+' '+str(((1+marker2_index)*10)+9)+' black\n')
										elif valeur <= VVV8:
											OUT.write('abscisse '+nom1+' '+pos_mark1d+' '+pos_mark1f+' ordonnee '+rename(DICO[name_marker2],LISTE)+' '+str((1+marker2_index)*10)+' '+str(((1+marker2_index)*10)+9)+' red\n')
										elif valeur <= VVV7:
											OUT.write('abscisse '+nom1+' '+pos_mark1d+' '+pos_mark1f+' ordonnee '+rename(DICO[name_marker2],LISTE)+' '+str((1+marker2_index)*10)+' '+str(((1+marker2_index)*10)+9)+' pink\n')
										elif valeur <= VVV6:
											OUT.write('abscisse '+nom1+' '+pos_mark1d+' '+pos_mark1f+' ordonnee '+rename(DICO[name_marker2],LISTE)+' '+str((1+marker2_index)*10)+' '+str(((1+marker2_index)*10)+9)+' orange\n')
										elif valeur <= VVV5:
											OUT.write('abscisse '+nom1+' '+pos_mark1d+' '+pos_mark1f+' ordonnee '+rename(DICO[name_marker2],LISTE)+' '+str((1+marker2_index)*10)+' '+str(((1+marker2_index)*10)+9)+' yellow\n')
										elif valeur <= VVV4:
											OUT.write('abscisse '+nom1+' '+pos_mark1d+' '+pos_mark1f+' ordonnee '+rename(DICO[name_marker2],LISTE)+' '+str((1+marker2_index)*10)+' '+str(((1+marker2_index)*10)+9)+' lightgreen\n')
										elif valeur <= VVV3:
											OUT.write('abscisse '+nom1+' '+pos_mark1d+' '+pos_mark1f+' ordonnee '+rename(DICO[name_marker2],LISTE)+' '+str((1+marker2_index)*10)+' '+str(((1+marker2_index)*10)+9)+' green\n')
										elif valeur <= VVV2:
											OUT.write('abscisse '+nom1+' '+pos_mark1d+' '+pos_mark1f+' ordonnee '+rename(DICO[name_marker2],LISTE)+' '+str((1+marker2_index)*10)+' '+str(((1+marker2_index)*10)+9)+' lightskyblue\n')
										else :
											OUT.write('abscisse '+nom1+' '+pos_mark1d+' '+pos_mark1f+' ordonnee '+rename(DICO[name_marker2],LISTE)+' '+str((1+marker2_index)*10)+' '+str(((1+marker2_index)*10)+9)+' blue3\n')
								if TYPE == 'COR':
									if valeur <= VV9:
										if valeur <= VV1:
											OUT.write('ordonnee '+nom1+' '+pos_mark1d+' '+pos_mark1f+' abscisse '+rename(DICO[name_marker2],LISTE)+' '+str((1+marker2_index)*10)+' '+str(((1+marker2_index)*10)+9)+' black\n')
										elif valeur <= VV2:
											OUT.write('ordonnee '+nom1+' '+pos_mark1d+' '+pos_mark1f+' abscisse '+rename(DICO[name_marker2],LISTE)+' '+str((1+marker2_index)*10)+' '+str(((1+marker2_index)*10)+9)+' red\n')
										elif valeur <= VV3:
											OUT.write('ordonnee '+nom1+' '+pos_mark1d+' '+pos_mark1f+' abscisse '+rename(DICO[name_marker2],LISTE)+' '+str((1+marker2_index)*10)+' '+str(((1+marker2_index)*10)+9)+' pink\n')
										elif valeur <= VV4:
											OUT.write('ordonnee '+nom1+' '+pos_mark1d+' '+pos_mark1f+' abscisse '+rename(DICO[name_marker2],LISTE)+' '+str((1+marker2_index)*10)+' '+str(((1+marker2_index)*10)+9)+' orange\n')
										elif valeur <= VV5:
											OUT.write('ordonnee '+nom1+' '+pos_mark1d+' '+pos_mark1f+' abscisse '+rename(DICO[name_marker2],LISTE)+' '+str((1+marker2_index)*10)+' '+str(((1+marker2_index)*10)+9)+' yellow\n')
										elif valeur <= VV6:
											OUT.write('ordonnee '+nom1+' '+pos_mark1d+' '+pos_mark1f+' abscisse '+rename(DICO[name_marker2],LISTE)+' '+str((1+marker2_index)*10)+' '+str(((1+marker2_index)*10)+9)+' lightgreen\n')
										elif valeur <= VV7:
											OUT.write('ordonnee '+nom1+' '+pos_mark1d+' '+pos_mark1f+' abscisse '+rename(DICO[name_marker2],LISTE)+' '+str((1+marker2_index)*10)+' '+str(((1+marker2_index)*10)+9)+' green\n')
										elif valeur <= VV8:
											OUT.write('ordonnee '+nom1+' '+pos_mark1d+' '+pos_mark1f+' abscisse '+rename(DICO[name_marker2],LISTE)+' '+str((1+marker2_index)*10)+' '+str(((1+marker2_index)*10)+9)+' lightskyblue\n')
										else :
											OUT.write('ordonnee '+nom1+' '+pos_mark1d+' '+pos_mark1f+' abscisse '+rename(DICO[name_marker2],LISTE)+' '+str((1+marker2_index)*10)+' '+str(((1+marker2_index)*10)+9)+' blue3\n')
							else:#region colonne avant region ligne
								if TYPE == 'LOD' or TYPE == 'COR':
									if valeur >= V9:
										if valeur >= V1:
											OUT.write('abscisse '+rename(DICO[name_marker2],LISTE)+' '+str((1+marker2_index)*10)+' '+str(((1+marker2_index)*10)+9)+' ordonnee '+nom1+' '+pos_mark1d+' '+pos_mark1f+' black\n')
										elif valeur >= V2:
											OUT.write('abscisse '+rename(DICO[name_marker2],LISTE)+' '+str((1+marker2_index)*10)+' '+str(((1+marker2_index)*10)+9)+' ordonnee '+nom1+' '+pos_mark1d+' '+pos_mark1f+' red\n')
										elif valeur >= V3:
											OUT.write('abscisse '+rename(DICO[name_marker2],LISTE)+' '+str((1+marker2_index)*10)+' '+str(((1+marker2_index)*10)+9)+' ordonnee '+nom1+' '+pos_mark1d+' '+pos_mark1f+' pink\n')
										elif valeur >= V4:
											OUT.write('abscisse '+rename(DICO[name_marker2],LISTE)+' '+str((1+marker2_index)*10)+' '+str(((1+marker2_index)*10)+9)+' ordonnee '+nom1+' '+pos_mark1d+' '+pos_mark1f+' orange\n')
										elif valeur >= V5:
											OUT.write('abscisse '+rename(DICO[name_marker2],LISTE)+' '+str((1+marker2_index)*10)+' '+str(((1+marker2_index)*10)+9)+' ordonnee '+nom1+' '+pos_mark1d+' '+pos_mark1f+' yellow\n')
										elif valeur >= V6:
											OUT.write('abscisse '+rename(DICO[name_marker2],LISTE)+' '+str((1+marker2_index)*10)+' '+str(((1+marker2_index)*10)+9)+' ordonnee '+nom1+' '+pos_mark1d+' '+pos_mark1f+' lightgreen\n')
										elif valeur >= V7:
											OUT.write('abscisse '+rename(DICO[name_marker2],LISTE)+' '+str((1+marker2_index)*10)+' '+str(((1+marker2_index)*10)+9)+' ordonnee '+nom1+' '+pos_mark1d+' '+pos_mark1f+' green\n')
										elif valeur >= V8:
											OUT.write('abscisse '+rename(DICO[name_marker2],LISTE)+' '+str((1+marker2_index)*10)+' '+str(((1+marker2_index)*10)+9)+' ordonnee '+nom1+' '+pos_mark1d+' '+pos_mark1f+' lightskyblue\n')
										else :
											OUT.write('abscisse '+rename(DICO[name_marker2],LISTE)+' '+str((1+marker2_index)*10)+' '+str(((1+marker2_index)*10)+9)+' ordonnee '+nom1+' '+pos_mark1d+' '+pos_mark1f+' blue3\n')
								elif TYPE == 'REC':
									if valeur <= VVV1:
										if valeur <= VVV9:
											OUT.write('abscisse '+rename(DICO[name_marker2],LISTE)+' '+str((1+marker2_index)*10)+' '+str(((1+marker2_index)*10)+9)+' ordonnee '+nom1+' '+pos_mark1d+' '+pos_mark1f+' black\n')
										elif valeur <= VVV8:
											OUT.write('abscisse '+rename(DICO[name_marker2],LISTE)+' '+str((1+marker2_index)*10)+' '+str(((1+marker2_index)*10)+9)+' ordonnee '+nom1+' '+pos_mark1d+' '+pos_mark1f+' red\n')
										elif valeur <= VVV7:
											OUT.write('abscisse '+rename(DICO[name_marker2],LISTE)+' '+str((1+marker2_index)*10)+' '+str(((1+marker2_index)*10)+9)+' ordonnee '+nom1+' '+pos_mark1d+' '+pos_mark1f+' pink\n')
										elif valeur <= VVV6:
											OUT.write('abscisse '+rename(DICO[name_marker2],LISTE)+' '+str((1+marker2_index)*10)+' '+str(((1+marker2_index)*10)+9)+' ordonnee '+nom1+' '+pos_mark1d+' '+pos_mark1f+' orange\n')
										elif valeur <= VVV5:
											OUT.write('abscisse '+rename(DICO[name_marker2],LISTE)+' '+str((1+marker2_index)*10)+' '+str(((1+marker2_index)*10)+9)+' ordonnee '+nom1+' '+pos_mark1d+' '+pos_mark1f+' yellow\n')
										elif valeur <= VVV4:
											OUT.write('abscisse '+rename(DICO[name_marker2],LISTE)+' '+str((1+marker2_index)*10)+' '+str(((1+marker2_index)*10)+9)+' ordonnee '+nom1+' '+pos_mark1d+' '+pos_mark1f+' lightgreen\n')
										elif valeur <= VVV3:
											OUT.write('abscisse '+rename(DICO[name_marker2],LISTE)+' '+str((1+marker2_index)*10)+' '+str(((1+marker2_index)*10)+9)+' ordonnee '+nom1+' '+pos_mark1d+' '+pos_mark1f+' green\n')
										elif valeur <= VVV2:
											OUT.write('abscisse '+rename(DICO[name_marker2],LISTE)+' '+str((1+marker2_index)*10)+' '+str(((1+marker2_index)*10)+9)+' ordonnee '+nom1+' '+pos_mark1d+' '+pos_mark1f+' lightskyblue\n')
										else :
											OUT.write('abscisse '+rename(DICO[name_marker2],LISTE)+' '+str((1+marker2_index)*10)+' '+str(((1+marker2_index)*10)+9)+' ordonnee '+nom1+' '+pos_mark1d+' '+pos_mark1f+' blue3\n')
								if TYPE == 'COR':
									if valeur <= VV9:
										if valeur <= VV1:
											OUT.write('ordonnee '+rename(DICO[name_marker2],LISTE)+' '+str((1+marker2_index)*10)+' '+str(((1+marker2_index)*10)+9)+' abscisse '+nom1+' '+pos_mark1d+' '+pos_mark1f+' black\n')
										elif valeur <= VV2:
											OUT.write('ordonnee '+rename(DICO[name_marker2],LISTE)+' '+str((1+marker2_index)*10)+' '+str(((1+marker2_index)*10)+9)+' abscisse '+nom1+' '+pos_mark1d+' '+pos_mark1f+' red\n')
										elif valeur <= VV3:
											OUT.write('ordonnee '+rename(DICO[name_marker2],LISTE)+' '+str((1+marker2_index)*10)+' '+str(((1+marker2_index)*10)+9)+' abscisse '+nom1+' '+pos_mark1d+' '+pos_mark1f+' pink\n')
										elif valeur <= VV4:
											OUT.write('ordonnee '+rename(DICO[name_marker2],LISTE)+' '+str((1+marker2_index)*10)+' '+str(((1+marker2_index)*10)+9)+' abscisse '+nom1+' '+pos_mark1d+' '+pos_mark1f+' orange\n')
										elif valeur <= VV5:
											OUT.write('ordonnee '+rename(DICO[name_marker2],LISTE)+' '+str((1+marker2_index)*10)+' '+str(((1+marker2_index)*10)+9)+' abscisse '+nom1+' '+pos_mark1d+' '+pos_mark1f+' yellow\n')
										elif valeur <= VV6:
											OUT.write('ordonnee '+rename(DICO[name_marker2],LISTE)+' '+str((1+marker2_index)*10)+' '+str(((1+marker2_index)*10)+9)+' abscisse '+nom1+' '+pos_mark1d+' '+pos_mark1f+' lightgreen\n')
										elif valeur <= VV7:
											OUT.write('ordonnee '+rename(DICO[name_marker2],LISTE)+' '+str((1+marker2_index)*10)+' '+str(((1+marker2_index)*10)+9)+' abscisse '+nom1+' '+pos_mark1d+' '+pos_mark1f+' green\n')
										elif valeur <= VV8:
											OUT.write('ordonnee '+rename(DICO[name_marker2],LISTE)+' '+str((1+marker2_index)*10)+' '+str(((1+marker2_index)*10)+9)+' abscisse '+nom1+' '+pos_mark1d+' '+pos_mark1f+' lightskyblue\n')
										else :
											OUT.write('ordonnee '+rename(DICO[name_marker2],LISTE)+' '+str((1+marker2_index)*10)+' '+str(((1+marker2_index)*10)+9)+' abscisse '+nom1+' '+pos_mark1d+' '+pos_mark1f+' blue3\n')
			i += 1
	OUT.flush()

def __main__():
	#Parse Command Line
	parser = optparse.OptionParser(usage="python %prog [options]\n\nProgram designed by Guillaume MARTIN : guillaume.martin@cirad.fr\n\n"
	"This program takes a matrix containing pairwise statistics between markers and a file containing markers order and plot pairwise marker statistics in a dot-plot like picture.\n")
	# Wrapper options. 
	parser.add_option( '', '--mat', dest='mat', default=None, help='Matrix file containing the pairwise information.')
	parser.add_option( '', '--order', dest='order', default=None, help='Table file containing marker ordered. column 1 : marker name, column 2 : identifier (ex: chr number). Markers should be grouped by identifier and ordered based on their position in the group. The file can contain additional columns (example position on identifier) that will not be used.')
	parser.add_option( '', '--type', dest='type', default='LOD', help='Type of statistics in the matrix : LOD (odd ratio), REC (recombination), COR (correlation)')
	parser.add_option( '', '--png', dest='png', default='dot-plot.png', help='Output file name. [default: %default]')
	(options, args) = parser.parse_args()
	
	if options.order == None:
		sys.exit('--order argument is missing')
	if options.mat == None:
		sys.exit('--mat argument is missing')
		
	ScriptPath = os.path.dirname(sys.argv[0])
	
	loca_programs = ConfigParser.RawConfigParser()
	loca_programs.read(ScriptPath+'/loca_programs.conf')
	
	##
	os.system('echo "Loading marker order"')
	file = open(options.order)
	dico = {}
	dico_inv = {}
	liste_chr = []
	liste_mark = []
	for line in file:
		data = line.split()
		if data != []:
			if data[1] in dico:
				dico[data[1]].append(data[0])
			else:
				liste_chr.append(data[1])
				dico[data[1]] = []
				dico[data[1]].append(data[0])
			liste_mark.append(data[0])
			dico_inv[data[0]] = data[1]
	file.close()
	
	# os.system('echo "Recording maximum and minimum"')
	min_max_val = min_max(options.mat)
	
	
	# os.system('echo "Creating the file for orthodotter"')
	temp = tempfile.NamedTemporaryFile()
	format_ortho(options.mat, dico_inv, liste_chr, liste_mark, min_max_val[0], min_max_val[1], options.type, temp)
	temp.flush()
	
	# os.system('echo "running orthodotter.pl"')
	qs=os.popen('wc -l '+options.order)
	value = ''
	for n in qs:
		value = int(n.split()[0])
	value = (value*3) + 100
	
	ortho = '%s -f %s  -toPlot ordonnee:abscisse -x %s -y %s -bg white -o %s -r 3 -fSize 2' % (ScriptPath+loca_programs.get('Programs','orthodotter'), temp.name, str(value), str(value), options.png)
	run_job(ortho, 'Error when running orthodotter:')
	temp.close()
	
if __name__ == "__main__": __main__()