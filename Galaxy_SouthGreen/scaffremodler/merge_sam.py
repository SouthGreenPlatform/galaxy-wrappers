
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

def CIGAR(mot, debut):
	if 'N' in mot:
		sys.exit('N in mot')
	elif 'P' in mot:
		sys.exit('P in mot')
	elif '=' in mot:
		sys.exit('= in mot')
	elif 'X' in mot:
		sys.exit('X in mot')
	else:
		mot_split = mot.replace('M',' M ').replace('D',' D ').replace('I',' I ').replace('S',' S ').replace('H',' H ').split()
	total = 0
	while mot_split:
		if mot_split[1] == 'M' or mot_split[1] == 'D':
			total = total + int(mot_split[0])
		del mot_split[0]
		del mot_split[0]
	return debut + (total-1)

def find_info(LINE, YT):
	mot = []
	for n in LINE:
		liste = n.split(':')
		if liste[0] == 'AS':
			mot.append(n)
		elif liste[0] == 'XS':
			mot.append(n)
	mot.append(YT)
	return mot


def merge2sam(FILE1, FILE2, OUT, MIN, MAX, OR):
	mapped_pair = 0
	mapped_single = 0
	unmapped = 0
	file1 = open(FILE1)
	file2 = open(FILE2)
	boucle = 1
	outfile = open(OUT, 'w')
	line1 = file1.readline()
	line2 = file2.readline()
	while line1:
		data1 = line1.split()
		data2 = line2.split()
		if data1[0][0] == '@':
			outfile.write(line1)
		else:
			while data1[1] == '256' or data1[1] == '272':
				line1 = file1.readline()
				data1 = line1.split()
			while data2[1] == '256' or data2[1] == '272':
				line2 = file2.readline()
				data2 = line2.split()
			if data1[0].replace('/1','').replace('/2','') != data2[0].replace('/1','').replace('/2','') and data1[0].replace('_1','').replace('_2','') != data2[0].replace('_1','').replace('_2',''):
				mot = 'Probleme in the mapping : '+data1[0]+' and '+data2[0]+' are different. Read mates should be identified with /1 and /2'
				sys.exit(mot)
			else:
				if data1[1] == '4' and data2[1] == '4': #reads unmapped
					unmapped += 1
					outfile.write(data1[0].replace('/1','').replace('/2','').replace('_1','').replace('_2','')+'\t77\t'+data1[2]+'\t'+data1[3]+'\t'+data1[4]+'\t'+data1[5]+'\t'+data1[6]+'\t'+data1[7]+'\t'+data1[8]+'\t'+data1[9]+'\t'+data1[10]+'\t'+'\t'.join(find_info(data1[11:],'YT:Z:UP'))+'\n')
					outfile.write(data2[0].replace('/1','').replace('/2','').replace('_1','').replace('_2','')+'\t141\t'+data2[2]+'\t'+data2[3]+'\t'+data2[4]+'\t'+data2[5]+'\t'+data2[6]+'\t'+data2[7]+'\t'+data2[8]+'\t'+data2[9]+'\t'+data2[10]+'\t'+'\t'.join(find_info(data2[11:],'YT:Z:UP'))+'\n')
				elif data1[1] == '0' and data2[1] == '4': #mate 1 mappant en F
					mapped_single += 1
					outfile.write(data1[0].replace('/1','').replace('/2','').replace('_1','').replace('_2','')+'\t73\t'+data1[2]+'\t'+data1[3]+'\t'+data1[4]+'\t'+data1[5]+'\t=\t'+data1[3]+'\t0\t'+data1[9]+'\t'+data1[10]+'\t'+'\t'.join(find_info(data1[11:],'YT:Z:UP'))+'\n')
					outfile.write(data2[0].replace('/1','').replace('/2','').replace('_1','').replace('_2','')+'\t133\t'+data1[2]+'\t'+data1[3]+'\t'+data2[4]+'\t'+data2[5]+'\t=\t'+data1[3]+'\t0\t'+data2[9]+'\t'+data2[10]+'\t'+'\t'.join(find_info(data2[11:],'YT:Z:UP'))+'\n')
				elif data1[1] == '16' and data2[1] == '4': #mate 1 mappant en R
					mapped_single += 1
					outfile.write(data1[0].replace('/1','').replace('/2','').replace('_1','').replace('_2','')+'\t89\t'+data1[2]+'\t'+data1[3]+'\t'+data1[4]+'\t'+data1[5]+'\t=\t'+data1[3]+'\t0\t'+data1[9]+'\t'+data1[10]+'\t'+'\t'.join(find_info(data1[11:],'YT:Z:UP'))+'\n')
					outfile.write(data2[0].replace('/1','').replace('/2','').replace('_1','').replace('_2','')+'\t165\t'+data1[2]+'\t'+data1[3]+'\t'+data2[4]+'\t'+data2[5]+'\t=\t'+data1[3]+'\t0\t'+data2[9]+'\t'+data2[10]+'\t'+'\t'.join(find_info(data2[11:],'YT:Z:UP'))+'\n')
				elif data1[1] == '4' and data2[1] == '0':#mate 2 mappant en F
					mapped_single += 1
					outfile.write(data2[0].replace('/1','').replace('/2','').replace('_1','').replace('_2','')+'\t137\t'+data2[2]+'\t'+data2[3]+'\t'+data2[4]+'\t'+data2[5]+'\t=\t'+data2[3]+'\t0\t'+data2[9]+'\t'+data2[10]+'\t'+'\t'.join(find_info(data2[11:],'YT:Z:UP'))+'\n')
					outfile.write(data1[0].replace('/1','').replace('/2','').replace('_1','').replace('_2','')+'\t69\t'+data2[2]+'\t'+data2[3]+'\t'+data1[4]+'\t'+data1[5]+'\t=\t'+data2[3]+'\t0\t'+data1[9]+'\t'+data1[10]+'\t'+'\t'.join(find_info(data1[11:],'YT:Z:UP'))+'\n')
				elif data1[1] == '4' and data2[1] == '16':#mate 2 mappant en R
					mapped_single += 1
					outfile.write(data2[0].replace('/1','').replace('/2','').replace('_1','').replace('_2','')+'\t153\t'+data2[2]+'\t'+data2[3]+'\t'+data2[4]+'\t'+data2[5]+'\t=\t'+data2[3]+'\t0\t'+data2[9]+'\t'+data2[10]+'\t'+'\t'.join(find_info(data2[11:],'YT:Z:UP'))+'\n')
					outfile.write(data1[0].replace('/1','').replace('/2','').replace('_1','').replace('_2','')+'\t101\t'+data2[2]+'\t'+data2[3]+'\t'+data1[4]+'\t'+data1[5]+'\t=\t'+data2[3]+'\t0\t'+data1[9]+'\t'+data1[10]+'\t'+'\t'.join(find_info(data1[11:],'YT:Z:UP'))+'\n')
				elif data1[1] == '0' and data2[1] == '0':#mate FF
					mapped_pair += 1
					if data1[2] != data2[2]:#sur des chromosomes differents
						outfile.write(data1[0].replace('/1','').replace('/2','').replace('_1','').replace('_2','')+'\t65\t'+data1[2]+'\t'+data1[3]+'\t'+data1[4]+'\t'+data1[5]+'\t'+data2[2]+'\t'+data2[3]+'\t0\t'+data1[9]+'\t'+data1[10]+'\t'+'\t'.join(find_info(data1[11:],'YT:Z:DP'))+'\n')
						outfile.write(data2[0].replace('/1','').replace('/2','').replace('_1','').replace('_2','')+'\t129\t'+data2[2]+'\t'+data2[3]+'\t'+data2[4]+'\t'+data2[5]+'\t'+data1[2]+'\t'+data1[3]+'\t0\t'+data2[9]+'\t'+data2[10]+'\t'+'\t'.join(find_info(data2[11:],'YT:Z:DP'))+'\n')
					else:
						if int(data1[3]) < int(data2[3]):
							outfile.write(data1[0].replace('/1','').replace('/2','').replace('_1','').replace('_2','')+'\t65\t'+data1[2]+'\t'+data1[3]+'\t'+data1[4]+'\t'+data1[5]+'\t=\t'+data2[3]+'\t'+str((CIGAR(data2[5],int(data2[3]))-int(data1[3]))+1)+'\t'+data1[9]+'\t'+data1[10]+'\t'+'\t'.join(find_info(data1[11:],'YT:Z:DP'))+'\n')
							outfile.write(data2[0].replace('/1','').replace('/2','').replace('_1','').replace('_2','')+'\t129\t'+data2[2]+'\t'+data2[3]+'\t'+data2[4]+'\t'+data2[5]+'\t=\t'+data1[3]+'\t-'+str((CIGAR(data2[5],int(data2[3]))-int(data1[3]))+1)+'\t'+data2[9]+'\t'+data2[10]+'\t'+'\t'.join(find_info(data2[11:],'YT:Z:DP'))+'\n')
						else:
							outfile.write(data1[0].replace('/1','').replace('/2','').replace('_1','').replace('_2','')+'\t65\t'+data1[2]+'\t'+data1[3]+'\t'+data1[4]+'\t'+data1[5]+'\t=\t'+data2[3]+'\t-'+str((CIGAR(data1[5],int(data1[3]))-int(data2[3]))+1)+'\t'+data1[9]+'\t'+data1[10]+'\t'+'\t'.join(find_info(data1[11:],'YT:Z:DP'))+'\n')
							outfile.write(data2[0].replace('/1','').replace('/2','').replace('_1','').replace('_2','')+'\t129\t'+data2[2]+'\t'+data2[3]+'\t'+data2[4]+'\t'+data2[5]+'\t=\t'+data1[3]+'\t'+str((CIGAR(data1[5],int(data1[3]))-int(data2[3]))+1)+'\t'+data2[9]+'\t'+data2[10]+'\t'+'\t'.join(find_info(data2[11:],'YT:Z:DP'))+'\n')
				elif data1[1] == '16' and data2[1] == '16':#mate RR
					mapped_pair += 1
					if data1[2] != data2[2]:#sur des chromosomes differents
						outfile.write(data1[0].replace('/1','').replace('/2','').replace('_1','').replace('_2','')+'\t113\t'+data1[2]+'\t'+data1[3]+'\t'+data1[4]+'\t'+data1[5]+'\t'+data2[2]+'\t'+data2[3]+'\t0\t'+data1[9]+'\t'+data1[10]+'\t'+'\t'.join(find_info(data1[11:],'YT:Z:DP'))+'\n')
						outfile.write(data2[0].replace('/1','').replace('/2','').replace('_1','').replace('_2','')+'\t177\t'+data2[2]+'\t'+data2[3]+'\t'+data2[4]+'\t'+data2[5]+'\t'+data1[2]+'\t'+data1[3]+'\t0\t'+data2[9]+'\t'+data2[10]+'\t'+'\t'.join(find_info(data2[11:],'YT:Z:DP'))+'\n')
					else:
						if int(data1[3]) < int(data2[3]):
							outfile.write(data1[0].replace('/1','').replace('/2','').replace('_1','').replace('_2','')+'\t113\t'+data1[2]+'\t'+data1[3]+'\t'+data1[4]+'\t'+data1[5]+'\t=\t'+data2[3]+'\t'+str((CIGAR(data2[5],int(data2[3]))-int(data1[3]))+1)+'\t'+data1[9]+'\t'+data1[10]+'\t'+'\t'.join(find_info(data1[11:],'YT:Z:DP'))+'\n')
							outfile.write(data2[0].replace('/1','').replace('/2','').replace('_1','').replace('_2','')+'\t177\t'+data2[2]+'\t'+data2[3]+'\t'+data2[4]+'\t'+data2[5]+'\t=\t'+data1[3]+'\t-'+str((CIGAR(data2[5],int(data2[3]))-int(data1[3]))+1)+'\t'+data2[9]+'\t'+data2[10]+'\t'+'\t'.join(find_info(data2[11:],'YT:Z:DP'))+'\n')
						else:
							outfile.write(data1[0].replace('/1','').replace('/2','').replace('_1','').replace('_2','')+'\t113\t'+data1[2]+'\t'+data1[3]+'\t'+data1[4]+'\t'+data1[5]+'\t=\t'+data2[3]+'\t-'+str((CIGAR(data1[5],int(data1[3]))-int(data2[3]))+1)+'\t'+data1[9]+'\t'+data1[10]+'\t'+'\t'.join(find_info(data1[11:],'YT:Z:DP'))+'\n')
							outfile.write(data2[0].replace('/1','').replace('/2','').replace('_1','').replace('_2','')+'\t177\t'+data2[2]+'\t'+data2[3]+'\t'+data2[4]+'\t'+data2[5]+'\t=\t'+data1[3]+'\t'+str((CIGAR(data1[5],int(data1[3]))-int(data2[3]))+1)+'\t'+data2[9]+'\t'+data2[10]+'\t'+'\t'.join(find_info(data2[11:],'YT:Z:DP'))+'\n')
				elif data1[1] == '16' and data2[1] == '0':#mate1 R, mate2 F
					mapped_pair += 1
					if data1[2] != data2[2]:#sur des chromosomes differents
						outfile.write(data1[0].replace('/1','').replace('/2','').replace('_1','').replace('_2','')+'\t81\t'+data1[2]+'\t'+data1[3]+'\t'+data1[4]+'\t'+data1[5]+'\t'+data2[2]+'\t'+data2[3]+'\t0\t'+data1[9]+'\t'+data1[10]+'\t'+'\t'.join(find_info(data1[11:],'YT:Z:DP'))+'\n')
						outfile.write(data2[0].replace('/1','').replace('/2','').replace('_1','').replace('_2','')+'\t161\t'+data2[2]+'\t'+data2[3]+'\t'+data2[4]+'\t'+data2[5]+'\t'+data1[2]+'\t'+data1[3]+'\t0\t'+data2[9]+'\t'+data2[10]+'\t'+'\t'.join(find_info(data2[11:],'YT:Z:DP'))+'\n')
					elif OR == 'rf':#pour une bonne orientation RF
						if int(data1[3]) < int(data2[3]) and ((CIGAR(data2[5],int(data2[3]))-int(data1[3]))+1) >= MIN and ((CIGAR(data2[5],int(data2[3]))-int(data1[3]))+1) <= MAX:#well mapped
							outfile.write(data1[0].replace('/1','').replace('/2','').replace('_1','').replace('_2','')+'\t83\t'+data1[2]+'\t'+data1[3]+'\t'+data1[4]+'\t'+data1[5]+'\t=\t'+data2[3]+'\t'+str((CIGAR(data2[5],int(data2[3]))-int(data1[3]))+1)+'\t'+data1[9]+'\t'+data1[10]+'\t'+'\t'.join(find_info(data1[11:],'YT:Z:CP'))+'\n')
							outfile.write(data2[0].replace('/1','').replace('/2','').replace('_1','').replace('_2','')+'\t163\t'+data2[2]+'\t'+data2[3]+'\t'+data2[4]+'\t'+data2[5]+'\t=\t'+data1[3]+'\t-'+str((CIGAR(data2[5],int(data2[3]))-int(data1[3]))+1)+'\t'+data2[9]+'\t'+data2[10]+'\t'+'\t'.join(find_info(data2[11:],'YT:Z:CP'))+'\n')
						else:
							if int(data1[3]) < int(data2[3]):#wrong insert size
								outfile.write(data1[0].replace('/1','').replace('/2','').replace('_1','').replace('_2','')+'\t81\t'+data1[2]+'\t'+data1[3]+'\t'+data1[4]+'\t'+data1[5]+'\t=\t'+data2[3]+'\t'+str((CIGAR(data2[5],int(data2[3]))-int(data1[3]))+1)+'\t'+data1[9]+'\t'+data1[10]+'\t'+'\t'.join(find_info(data1[11:],'YT:Z:DP'))+'\n')
								outfile.write(data2[0].replace('/1','').replace('/2','').replace('_1','').replace('_2','')+'\t161\t'+data2[2]+'\t'+data2[3]+'\t'+data2[4]+'\t'+data2[5]+'\t=\t'+data1[3]+'\t-'+str((CIGAR(data2[5],int(data2[3]))-int(data1[3]))+1)+'\t'+data2[9]+'\t'+data2[10]+'\t'+'\t'.join(find_info(data2[11:],'YT:Z:DP'))+'\n')
							else:#read FR
								outfile.write(data1[0].replace('/1','').replace('/2','').replace('_1','').replace('_2','')+'\t81\t'+data1[2]+'\t'+data1[3]+'\t'+data1[4]+'\t'+data1[5]+'\t=\t'+data2[3]+'\t-'+str((CIGAR(data1[5],int(data1[3]))-int(data2[3]))+1)+'\t'+data1[9]+'\t'+data1[10]+'\t'+'\t'.join(find_info(data1[11:],'YT:Z:DP'))+'\n')
								outfile.write(data2[0].replace('/1','').replace('/2','').replace('_1','').replace('_2','')+'\t161\t'+data2[2]+'\t'+data2[3]+'\t'+data2[4]+'\t'+data2[5]+'\t=\t'+data1[3]+'\t'+str((CIGAR(data1[5],int(data1[3]))-int(data2[3]))+1)+'\t'+data2[9]+'\t'+data2[10]+'\t'+'\t'.join(find_info(data2[11:],'YT:Z:DP'))+'\n')
					elif OR == 'fr':#pour une bonne orientation FR
						if int(data1[3]) > int(data2[3]) and ((CIGAR(data1[5],int(data1[3]))-int(data2[3]))+1) >= MIN and ((CIGAR(data1[5],int(data1[3]))-int(data2[3]))+1) <= MAX:
							outfile.write(data1[0].replace('/1','').replace('/2','').replace('_1','').replace('_2','')+'\t83\t'+data1[2]+'\t'+data1[3]+'\t'+data1[4]+'\t'+data1[5]+'\t=\t'+data2[3]+'\t-'+str((CIGAR(data1[5],int(data1[3]))-int(data2[3]))+1)+'\t'+data1[9]+'\t'+data1[10]+'\t'+'\t'.join(find_info(data1[11:],'YT:Z:CP'))+'\n')
							outfile.write(data2[0].replace('/1','').replace('/2','').replace('_1','').replace('_2','')+'\t163\t'+data2[2]+'\t'+data2[3]+'\t'+data2[4]+'\t'+data2[5]+'\t=\t'+data1[3]+'\t'+str((CIGAR(data1[5],int(data1[3]))-int(data2[3]))+1)+'\t'+data2[9]+'\t'+data2[10]+'\t'+'\t'.join(find_info(data2[11:],'YT:Z:CP'))+'\n')
						else:
							if int(data1[3]) < int(data2[3]):
								outfile.write(data1[0].replace('/1','').replace('/2','').replace('_1','').replace('_2','')+'\t81\t'+data1[2]+'\t'+data1[3]+'\t'+data1[4]+'\t'+data1[5]+'\t=\t'+data2[3]+'\t'+str((CIGAR(data2[5],int(data2[3]))-int(data1[3]))+1)+'\t'+data1[9]+'\t'+data1[10]+'\t'+'\t'.join(find_info(data1[11:],'YT:Z:DP'))+'\n')
								outfile.write(data2[0].replace('/1','').replace('/2','').replace('_1','').replace('_2','')+'\t161\t'+data2[2]+'\t'+data2[3]+'\t'+data2[4]+'\t'+data2[5]+'\t=\t'+data1[3]+'\t-'+str((CIGAR(data2[5],int(data2[3]))-int(data1[3]))+1)+'\t'+data2[9]+'\t'+data2[10]+'\t'+'\t'.join(find_info(data2[11:],'YT:Z:DP'))+'\n')
							else:
								outfile.write(data1[0].replace('/1','').replace('/2','').replace('_1','').replace('_2','')+'\t81\t'+data1[2]+'\t'+data1[3]+'\t'+data1[4]+'\t'+data1[5]+'\t=\t'+data2[3]+'\t-'+str((CIGAR(data1[5],int(data1[3]))-int(data2[3]))+1)+'\t'+data1[9]+'\t'+data1[10]+'\t'+'\t'.join(find_info(data1[11:],'YT:Z:DP'))+'\n')
								outfile.write(data2[0].replace('/1','').replace('/2','').replace('_1','').replace('_2','')+'\t161\t'+data2[2]+'\t'+data2[3]+'\t'+data2[4]+'\t'+data2[5]+'\t=\t'+data1[3]+'\t'+str((CIGAR(data1[5],int(data1[3]))-int(data2[3]))+1)+'\t'+data2[9]+'\t'+data2[10]+'\t'+'\t'.join(find_info(data2[11:],'YT:Z:DP'))+'\n')
					else:
						sys.exit('bug')
				elif data1[1] == '0' and data2[1] == '16':#mate1 F, mate2 R
					mapped_pair += 1
					if data1[2] != data2[2]:#sur des chromosomes differents
						outfile.write(data1[0].replace('/1','').replace('/2','').replace('_1','').replace('_2','')+'\t97\t'+data1[2]+'\t'+data1[3]+'\t'+data1[4]+'\t'+data1[5]+'\t'+data2[2]+'\t'+data2[3]+'\t0\t'+data1[9]+'\t'+data1[10]+'\t'+'\t'.join(find_info(data1[11:],'YT:Z:DP'))+'\n')
						outfile.write(data2[0].replace('/1','').replace('/2','').replace('_1','').replace('_2','')+'\t145\t'+data2[2]+'\t'+data2[3]+'\t'+data2[4]+'\t'+data2[5]+'\t'+data1[2]+'\t'+data1[3]+'\t0\t'+data2[9]+'\t'+data2[10]+'\t'+'\t'.join(find_info(data2[11:],'YT:Z:DP'))+'\n')
					elif OR == 'rf':#pour une bonne orientation RF
						if int(data1[3]) > int(data2[3]) and ((CIGAR(data1[5],int(data1[3]))-int(data2[3]))+1) >= MIN and ((CIGAR(data1[5],int(data1[3]))-int(data2[3]))+1) <= MAX:
							outfile.write(data1[0].replace('/1','').replace('/2','').replace('_1','').replace('_2','')+'\t99\t'+data1[2]+'\t'+data1[3]+'\t'+data1[4]+'\t'+data1[5]+'\t=\t'+data2[3]+'\t-'+str((CIGAR(data1[5],int(data1[3]))-int(data2[3]))+1)+'\t'+data1[9]+'\t'+data1[10]+'\t'+'\t'.join(find_info(data1[11:],'YT:Z:CP'))+'\n')
							outfile.write(data2[0].replace('/1','').replace('/2','').replace('_1','').replace('_2','')+'\t147\t'+data2[2]+'\t'+data2[3]+'\t'+data2[4]+'\t'+data2[5]+'\t=\t'+data1[3]+'\t'+str((CIGAR(data1[5],int(data1[3]))-int(data2[3]))+1)+'\t'+data2[9]+'\t'+data2[10]+'\t'+'\t'.join(find_info(data2[11:],'YT:Z:CP'))+'\n')
						else:
							if int(data1[3]) < int(data2[3]):
								outfile.write(data1[0].replace('/1','').replace('/2','').replace('_1','').replace('_2','')+'\t97\t'+data1[2]+'\t'+data1[3]+'\t'+data1[4]+'\t'+data1[5]+'\t=\t'+data2[3]+'\t'+str((CIGAR(data2[5],int(data2[3]))-int(data1[3]))+1)+'\t'+data1[9]+'\t'+data1[10]+'\t'+'\t'.join(find_info(data1[11:],'YT:Z:DP'))+'\n')
								outfile.write(data2[0].replace('/1','').replace('/2','').replace('_1','').replace('_2','')+'\t145\t'+data2[2]+'\t'+data2[3]+'\t'+data2[4]+'\t'+data2[5]+'\t=\t'+data1[3]+'\t-'+str((CIGAR(data2[5],int(data2[3]))-int(data1[3]))+1)+'\t'+data2[9]+'\t'+data2[10]+'\t'+'\t'.join(find_info(data2[11:],'YT:Z:DP'))+'\n')
							else:
								outfile.write(data1[0].replace('/1','').replace('/2','').replace('_1','').replace('_2','')+'\t97\t'+data1[2]+'\t'+data1[3]+'\t'+data1[4]+'\t'+data1[5]+'\t=\t'+data2[3]+'\t-'+str((CIGAR(data1[5],int(data1[3]))-int(data2[3]))+1)+'\t'+data1[9]+'\t'+data1[10]+'\t'+'\t'.join(find_info(data1[11:],'YT:Z:DP'))+'\n')
								outfile.write(data2[0].replace('/1','').replace('/2','').replace('_1','').replace('_2','')+'\t145\t'+data2[2]+'\t'+data2[3]+'\t'+data2[4]+'\t'+data2[5]+'\t=\t'+data1[3]+'\t'+str((CIGAR(data1[5],int(data1[3]))-int(data2[3]))+1)+'\t'+data2[9]+'\t'+data2[10]+'\t'+'\t'.join(find_info(data2[11:],'YT:Z:DP'))+'\n')
					elif OR == 'fr':#pour une bonne orientation FR
						if int(data1[3]) < int(data2[3]) and ((CIGAR(data2[5],int(data2[3]))-int(data1[3]))+1) >= MIN and ((CIGAR(data2[5],int(data2[3]))-int(data1[3]))+1) <= MAX:
							outfile.write(data1[0].replace('/1','').replace('/2','').replace('_1','').replace('_2','')+'\t99\t'+data1[2]+'\t'+data1[3]+'\t'+data1[4]+'\t'+data1[5]+'\t=\t'+data2[3]+'\t'+str((CIGAR(data2[5],int(data2[3]))-int(data1[3]))+1)+'\t'+data1[9]+'\t'+data1[10]+'\t'+'\t'.join(find_info(data1[11:],'YT:Z:CP'))+'\n')
							outfile.write(data2[0].replace('/1','').replace('/2','').replace('_1','').replace('_2','')+'\t147\t'+data2[2]+'\t'+data2[3]+'\t'+data2[4]+'\t'+data2[5]+'\t=\t'+data1[3]+'\t-'+str((CIGAR(data2[5],int(data2[3]))-int(data1[3]))+1)+'\t'+data2[9]+'\t'+data2[10]+'\t'+'\t'.join(find_info(data2[11:],'YT:Z:CP'))+'\n')
						else:
							if int(data1[3]) < int(data2[3]):
								outfile.write(data1[0].replace('/1','').replace('/2','').replace('_1','').replace('_2','')+'\t97\t'+data1[2]+'\t'+data1[3]+'\t'+data1[4]+'\t'+data1[5]+'\t=\t'+data2[3]+'\t'+str((CIGAR(data2[5],int(data2[3]))-int(data1[3]))+1)+'\t'+data1[9]+'\t'+data1[10]+'\t'+'\t'.join(find_info(data1[11:],'YT:Z:DP'))+'\n')
								outfile.write(data2[0].replace('/1','').replace('/2','').replace('_1','').replace('_2','')+'\t145\t'+data2[2]+'\t'+data2[3]+'\t'+data2[4]+'\t'+data2[5]+'\t=\t'+data1[3]+'\t-'+str((CIGAR(data2[5],int(data2[3]))-int(data1[3]))+1)+'\t'+data2[9]+'\t'+data2[10]+'\t'+'\t'.join(find_info(data2[11:],'YT:Z:DP'))+'\n')
							else:
								outfile.write(data1[0].replace('/1','').replace('/2','').replace('_1','').replace('_2','')+'\t97\t'+data1[2]+'\t'+data1[3]+'\t'+data1[4]+'\t'+data1[5]+'\t=\t'+data2[3]+'\t-'+str((CIGAR(data1[5],int(data1[3]))-int(data2[3]))+1)+'\t'+data1[9]+'\t'+data1[10]+'\t'+'\t'.join(find_info(data1[11:],'YT:Z:DP'))+'\n')
								outfile.write(data2[0].replace('/1','').replace('/2','').replace('_1','').replace('_2','')+'\t145\t'+data2[2]+'\t'+data2[3]+'\t'+data2[4]+'\t'+data2[5]+'\t=\t'+data1[3]+'\t'+str((CIGAR(data1[5],int(data1[3]))-int(data2[3]))+1)+'\t'+data2[9]+'\t'+data2[10]+'\t'+'\t'.join(find_info(data2[11:],'YT:Z:DP'))+'\n')
					else:
						sys.exit('bug')
				else:
					os.system('echo "'+data1[1]+' '+data2[1]+'"')
					sys.exit('Probleme in the formating of mapping file')
		line1 = file1.readline()
		line2 = file2.readline()
	outfile.close()
	os.system('echo "Mapped pair: '+str(mapped_pair)+'"')
	os.system('echo "Mapped single (mate1 or mate2): '+str(mapped_single)+'"')
	os.system('echo "Unmapped (mate1 and mate2): '+str(unmapped)+'"')
	


def __main__():
	#Parse Command Line
	parser = optparse.OptionParser(usage="python %prog [options]\n\nProgram designed by Guillaume MARTIN : guillaume.martin@cirad.fr")
	# Wrapper options.
	parser.add_option( '', '--file1', dest='file1', help='Mate1 sam file')
	parser.add_option( '', '--file2', dest='file2', help='Mate2 sam file')
	parser.add_option( '', '--out', dest='out', help='Output file name')
	parser.add_option( '', '--min', dest='min', help='minimal insert size to accept the pair as properly mapped')
	parser.add_option( '', '--max', dest='max', help='maximal insert size to accept the pair as properly mapped')
	parser.add_option( '', '--orient', dest='orient', help='Expected orientation of paired reads')
	(options, args) = parser.parse_args()
	
	
	merge2sam(options.file1, options.file2, options.out, int(options.min), int(options.max), options.orient)
	
if __name__ == "__main__": __main__()
