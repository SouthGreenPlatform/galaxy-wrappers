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

import optparse, os, shutil, subprocess, sys, tempfile, fileinput, ConfigParser, operator, time, math

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

def JMpwd2markid(PWD):
	dico = set()
	file = open(PWD)
	sur_id = 0
	for line in file:
		data = line.split()
		if data:
			if sur_id:
				dico.add(data[1])
			elif data[0] == 'locus' and data[1] == 'numbers:':
				sur_id = 1
	return dico

def JMpwd2matrix(PWD, LISTE_ID, OUT1, OUT2):
	#creating the table structure
	# os.system("echo '1)creating the table structure'")
	liste_vide = []
	for n in LISTE_ID:
		liste_vide.append('999999999')
	dic_LOD = {}
	dic_REC = {}
	for n in LISTE_ID:
		dic_LOD[n] = list(liste_vide)
		dic_REC[n] = list(liste_vide)
	
	#filling the table
	# os.system("echo '2)filling the table. It can takes time...'")
	file = open(PWD)
	maxlod = 0
	for line in file:
		data = line.split()
		if data != [] and not('name =' in line) and len(data) > 2:
			if data[0] != ';':
				if data[0] in LISTE_ID and  data[1] in LISTE_ID:
					dic_LOD[data[0]][LISTE_ID.index(data[1])] = float(data[3])
					dic_LOD[data[1]][LISTE_ID.index(data[0])] = float(data[3])
					dic_REC[data[0]][LISTE_ID.index(data[1])] = float(data[2])
					dic_REC[data[1]][LISTE_ID.index(data[0])] = float(data[2])
					if maxlod < float(data[3]):
						maxlod = float(data[3])
	file.close()
	
	# filling diag
	# os.system("echo '3)filling diag'")
	for n in LISTE_ID:
		dic_REC[n][LISTE_ID.index(n)] = 0
		dic_LOD[n][LISTE_ID.index(n)] = maxlod
			
	#writing results
	# os.system("echo '4)writing results'")
	OUT1.write('ID\t'+'\t'.join(LISTE_ID)+'\n')
	OUT2.write('ID\t'+'\t'.join(LISTE_ID)+'\n')
	for n in LISTE_ID:
		OUT1.write(n+'\t'+'\t'.join(map(str,dic_LOD[n]))+'\n')
		OUT2.write(n+'\t'+'\t'.join(map(str,dic_REC[n]))+'\n')
	OUT1.flush()
	OUT2.flush()

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
			min_val = min(val)
			max_val = max(val)
			# print min_val, max_val, data[0]
			if mini > min_val:
				mini = min_val
			if maxi < max_val:
				maxi = max_val
		# if i%1000 == 0:
			# os.system('echo "'+str(i)+'"')
	# os.system('echo "'+str(mini)+' '+str(maxi)+'"')
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
	VV1 = (MAX - (1*(MAX/40)))
	VV2 = (MAX - (2*(MAX/40)))
	VV3 = (MAX - (3*(MAX/40)))
	VV4 = (MAX - (4*(MAX/40)))
	VV5 = (MAX - (5*(MAX/40)))
	VV6 = (MAX - (6*(MAX/40)))
	VV7 = (MAX - (7*(MAX/40)))
	VV8 = (MAX - (8*(MAX/40)))
	VV9 = (MAX - (9*(MAX/40)))
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
									if valeur <= (MIN - (9*(MIN/10))):
										if valeur <= (MIN - (1*(MIN/10))):
											OUT.write('ordonnee '+nom1+' '+pos_mark1d+' '+pos_mark1f+' abscisse '+rename(DICO[name_marker2],LISTE)+' '+str((1+marker2_index)*10)+' '+str(((1+marker2_index)*10)+9)+' midnightblue\n')
										elif valeur <= (MIN - (2*(MIN/10))):
											OUT.write('ordonnee '+nom1+' '+pos_mark1d+' '+pos_mark1f+' abscisse '+rename(DICO[name_marker2],LISTE)+' '+str((1+marker2_index)*10)+' '+str(((1+marker2_index)*10)+9)+' navyblue\n')
										elif valeur <= (MIN - (3*(MIN/10))):
											OUT.write('ordonnee '+nom1+' '+pos_mark1d+' '+pos_mark1f+' abscisse '+rename(DICO[name_marker2],LISTE)+' '+str((1+marker2_index)*10)+' '+str(((1+marker2_index)*10)+9)+' darkblue\n')
										elif valeur <= (MIN - (4*(MIN/10))):
											OUT.write('ordonnee '+nom1+' '+pos_mark1d+' '+pos_mark1f+' abscisse '+rename(DICO[name_marker2],LISTE)+' '+str((1+marker2_index)*10)+' '+str(((1+marker2_index)*10)+9)+' blue4\n')
										elif valeur <= (MIN - (5*(MIN/10))):
											OUT.write('ordonnee '+nom1+' '+pos_mark1d+' '+pos_mark1f+' abscisse '+rename(DICO[name_marker2],LISTE)+' '+str((1+marker2_index)*10)+' '+str(((1+marker2_index)*10)+9)+' blue3\n')
										elif valeur <= (MIN - (6*(MIN/10))):
											OUT.write('ordonnee '+nom1+' '+pos_mark1d+' '+pos_mark1f+' abscisse '+rename(DICO[name_marker2],LISTE)+' '+str((1+marker2_index)*10)+' '+str(((1+marker2_index)*10)+9)+' blue2\n')
										elif valeur <= (MIN - (7*(MIN/10))):
											OUT.write('ordonnee '+nom1+' '+pos_mark1d+' '+pos_mark1f+' abscisse '+rename(DICO[name_marker2],LISTE)+' '+str((1+marker2_index)*10)+' '+str(((1+marker2_index)*10)+9)+' blue1\n')
										elif valeur <= (MIN - (8*(MIN/10))):
											OUT.write('ordonnee '+nom1+' '+pos_mark1d+' '+pos_mark1f+' abscisse '+rename(DICO[name_marker2],LISTE)+' '+str((1+marker2_index)*10)+' '+str(((1+marker2_index)*10)+9)+' lightskyblue\n')
										else :
											OUT.write('ordonnee '+nom1+' '+pos_mark1d+' '+pos_mark1f+' abscisse '+rename(DICO[name_marker2],LISTE)+' '+str((1+marker2_index)*10)+' '+str(((1+marker2_index)*10)+9)+' lightblue\n')
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
									if valeur <= (MIN - (9*(MIN/10))):
										if valeur <= (MIN - (1*(MIN/10))):
											OUT.write('ordonnee '+rename(DICO[name_marker2],LISTE)+' '+str((1+marker2_index)*10)+' '+str(((1+marker2_index)*10)+9)+' abscisse '+nom1+' '+pos_mark1d+' '+pos_mark1f+' midnightblue\n')
										elif valeur <= (MIN - (2*(MIN/10))):
											OUT.write('ordonnee '+rename(DICO[name_marker2],LISTE)+' '+str((1+marker2_index)*10)+' '+str(((1+marker2_index)*10)+9)+' abscisse '+nom1+' '+pos_mark1d+' '+pos_mark1f+' navyblue\n')
										elif valeur <= (MIN - (3*(MIN/10))):
											OUT.write('ordonnee '+rename(DICO[name_marker2],LISTE)+' '+str((1+marker2_index)*10)+' '+str(((1+marker2_index)*10)+9)+' abscisse '+nom1+' '+pos_mark1d+' '+pos_mark1f+' darkblue\n')
										elif valeur <= (MIN - (4*(MIN/10))):
											OUT.write('ordonnee '+rename(DICO[name_marker2],LISTE)+' '+str((1+marker2_index)*10)+' '+str(((1+marker2_index)*10)+9)+' abscisse '+nom1+' '+pos_mark1d+' '+pos_mark1f+' blue4\n')
										elif valeur <= (MIN - (5*(MIN/10))):
											OUT.write('ordonnee '+rename(DICO[name_marker2],LISTE)+' '+str((1+marker2_index)*10)+' '+str(((1+marker2_index)*10)+9)+' abscisse '+nom1+' '+pos_mark1d+' '+pos_mark1f+' blue3\n')
										elif valeur <= (MIN - (6*(MIN/10))):
											OUT.write('ordonnee '+rename(DICO[name_marker2],LISTE)+' '+str((1+marker2_index)*10)+' '+str(((1+marker2_index)*10)+9)+' abscisse '+nom1+' '+pos_mark1d+' '+pos_mark1f+' blue2\n')
										elif valeur <= (MIN - (7*(MIN/10))):
											OUT.write('ordonnee '+rename(DICO[name_marker2],LISTE)+' '+str((1+marker2_index)*10)+' '+str(((1+marker2_index)*10)+9)+' abscisse '+nom1+' '+pos_mark1d+' '+pos_mark1f+' blue1\n')
										elif valeur <= (MIN - (8*(MIN/10))):
											OUT.write('ordonnee '+rename(DICO[name_marker2],LISTE)+' '+str((1+marker2_index)*10)+' '+str(((1+marker2_index)*10)+9)+' abscisse '+nom1+' '+pos_mark1d+' '+pos_mark1f+' lightskyblue\n')
										else :
											OUT.write('ordonnee '+rename(DICO[name_marker2],LISTE)+' '+str((1+marker2_index)*10)+' '+str(((1+marker2_index)*10)+9)+' abscisse '+nom1+' '+pos_mark1d+' '+pos_mark1f+' lightblue\n')
			i += 1
	OUT.flush()

def matrix2ortho(MAT, ORDER, TYPE, OUT):
	##
	# os.system('echo "5)Creation of the file for orthodotter"')
	file = open(ORDER)
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
	
	##
	# os.system('echo "6)Recording maximum and minimum"')
	min_max_val = min_max(MAT)
	
	##
	# os.system('echo "Creating the file for orthodotter"')
	format_ortho(MAT, dico_inv, liste_chr, liste_mark, min_max_val[0], min_max_val[1], TYPE, OUT)

def rec2kosambi(MAT, OUT1, OUT2, ID):
	file = open(ID)
	dico = {}
	for line in file:
		data = line.split()
		if data:
			dico[data[0]] = data[1]
	file.close()
	outfile = open(OUT1,'w')
	outfile.write('@DARwin 5.0 - DIS\n')
	outfile2 = open(OUT2,'w')
	outfile2.write('@DARwin 5.0 - DON\t\n')
	file = open(MAT)
	line = file.readline()
	head = line.split()
	outfile.write(str(len(head[1:]))+'\n')
	outfile2.write(str(len(head[1:]))+'\t2\n')
	outfile2.write('unit\tname\tchrom\n')
	for n in head[1:-1]:
		outfile.write('\t'+str(head[1:].index(n)+1))
	for n in head[1:]:
		outfile2.write(str(head[1:].index(n)+1)+'\t'+n+'\t'+dico[n]+'\n')
	outfile.write('\n')
	outfile2.close()
	i = 0
	for line in file:
		data = line.split()
		if data:
			i += 1
			j = 0
			if i != 1:
				outfile.write(str(i))
				for n in data[1:]:
					j += 1
					if j < i:
						outfile.write('\t'+str( (math.log((1+ (2*float(n))) / (1- (2*float(n)))))/4 ))
				outfile.write('\n')
	outfile.close()

def __main__():
	#Parse Command Line
	parser = optparse.OptionParser(usage="python %prog [options]\n\nProgram designed by Guillaume MARTIN : guillaume.martin@cirad.fr\n\n"
	"This program take a pairwise file generated by JoinMap and return a dot-plot representing pairwise marker linkage of markers along "
	"chromosomes. Two additionnal file are generated that can be used with Darwin to construct tree.\n")
	# Wrapper options. 
	parser.add_option( '', '--pwd', dest='pwd', default=None, help='Pairwise file calculated by Joinmap.')
	parser.add_option( '', '--order', dest='order', default=None, help='Tabulated file containing in col1 : markers names and col2 : chromosome name. Markers should be grouped by identifier and ordered based on their position in the group. The file can contain additional columns (example position on identifier) that will not be used.')
	parser.add_option( '', '--type', dest='type', default='LOD', help='Type of pairwise statistics: LOD or REC. [default: %default]')
	parser.add_option( '', '--png', dest='png', default='dot-plot.png', help='Output name for the dot-plot (png file). [default: %default]')
	parser.add_option( '', '--dis', dest='dis', default='darwin.dis', help='Output name for the .dis file (for Darwin). [default: %default]')
	parser.add_option( '', '--don', dest='don', default='darwin.don', help='Output name for the .don file (for Darwin). [default: %default]')
	(options, args) = parser.parse_args()
	
	if options.order == None:
		sys.exit('--order argument is missing')
	if options.pwd == None:
		sys.exit('--pwd argument is missing')
	
	nom_pairwise = options.pwd
	
	ScriptPath = os.path.dirname(sys.argv[0])
	
	loca_programs = ConfigParser.RawConfigParser()
	loca_programs.read(ScriptPath+'/loca_programs.conf')
	
	
	qs=os.popen('wc -l '+options.order)
	value = ''
	for n in qs:
		value = int(n.split()[0])
	value = (value*3) + 100
	
	#creating a list of marker found in the pairwise file
	dico_id = JMpwd2markid(nom_pairwise)
	
	#A piece of verification
	file = open(options.order)
	liste_id = []
	for line in file:
		data = line.split()
		if data:
			if data[0] in liste_id:
				mot = 'The program exited without finishing: the marker '+data[0]+' found in more than once in '+options.order
				sys.exit(mot)
			liste_id.append(data[0])
			if not(data[0] in dico_id):
				mot = 'The program exited without finishing: the marker '+data[0]+' found in '+options.order+' is not in '+options.pwd
				sys.exit(mot)
	
	# creation of the matrix file
	temp2 = tempfile.NamedTemporaryFile()
	temp3 = tempfile.NamedTemporaryFile()
	JMpwd2matrix(nom_pairwise, liste_id, temp2, temp3)
	temp2.flush()
	temp3.flush()
	
	#creation of the file for orthodotter
	temp4 = tempfile.NamedTemporaryFile()
	if options.type == 'LOD':
		matrix2ortho (temp2.name,  options.order, options.type, temp4)
	elif options.type == 'REC':
		matrix2ortho (temp3.name,  options.order, options.type, temp4)
	else:
		mot = 'The program exited without finishing: the argument '+options.type+' is not recognized in --type'
		sys.exit(mot)
	temp4.flush()
	
	ortho = '%s -f %s  -toPlot ordonnee:abscisse -x %s -y %s -bg white -o %s -r 3 -fSize 2' % (ScriptPath+loca_programs.get('Programs','orthodotter'), temp4.name, str(value), str(value), options.png)
	run_job(ortho, 'Error when running orthodotter:')

	rec2kosambi(temp3.name, options.dis, options.don, options.order)
	
	temp2.close()
	temp3.close()
	temp4.close()
	
if __name__ == "__main__": __main__()