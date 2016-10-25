
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

def trie2discord_pair(LOCA_PROGRAMS, SAM, TYPE, SORT, MINI_DIS, MINI, MAXI, ORIENT, CHR, LISTE):
	if TYPE == 'bam':
		tmp = (tempfile.NamedTemporaryFile().name)+'.sam'
		if SORT == 'coordinate' or SORT == 'unsorted':
			bam2sam = '%s -jar %s SortSam INPUT=%s OUTPUT=%s SORT_ORDER=queryname QUIET=true MAX_RECORDS_IN_RAM=5000000 VERBOSITY=WARNING VALIDATION_STRINGENCY=SILENT' % (LOCA_PROGRAMS.get('Programs','java'), LOCA_PROGRAMS.get('Programs','picard-tool'), SAM, tmp)
		elif SORT == 'queryname':
			bam2sam = 'samtools view -h -o %s %s' % (tmp, SAM)
		else:
			mot = 'Unrecognized --sort option : '+SORT
			sys.exit(mot)
		run_job (bam2sam, 'Error in bam2sam conversion:')
		trielinebyline(tmp, MINI_DIS, MINI, MAXI, ORIENT, CHR, LISTE)
		# print tmp
		os.remove(tmp)
	else:
		if SORT == 'coordinate' or SORT == 'unsorted':
			tmp = (tempfile.NamedTemporaryFile().name)+'.sam'
			Sorting = '%s -jar %s SortSam INPUT=%s OUTPUT=%s SORT_ORDER=queryname QUIET=true MAX_RECORDS_IN_RAM=5000000 VERBOSITY=WARNING VALIDATION_STRINGENCY=SILENT' % (LOCA_PROGRAMS.get('Programs','java'), LOCA_PROGRAMS.get('Programs','picard-tool'), SAM, tmp)
			run_job (Sorting, 'Error in samtools :')
			trielinebyline(SAM, MINI_DIS, MINI, MAXI, ORIENT, CHR, LISTE)
			# print tmp
			os.remove(tmp)
		elif SORT == 'queryname':
			trielinebyline(SAM, MINI_DIS, MINI, MAXI, ORIENT, CHR, LISTE)
		else:
			mot = 'Unrecognized --sort option : '+SORT
			sys.exit(mot)

def trielinebyline(SAM, MINI_DIS, MINI, MAXI, ORIENT, CHR, LISTE):
	############################################
	#recording chromosome order
	############################################
	liste_chr = []
	file = open(CHR)
	for line in file:
		data = line.split()
		if data:
			liste_chr.append(data[0])
	file.close()

	############################################
	#Creating the liste of discordant reads
	############################################
	file = open(SAM)
	outfile = open(LISTE,'w')
	l1 = file.readline().split()
	while l1[0][0] == '@':
		l1 = file.readline().split()
	l2 = file.readline().split()
	while l1:
		if ORIENT == 'rf':
			outfile.write('\t'.join(trie_discord_rf(l1, l2, MINI_DIS, MINI, MAXI, liste_chr))+'\n')
		elif ORIENT == 'fr':
			outfile.write('\t'.join(trie_discord_fr(l1, l2, MINI_DIS, MINI, MAXI, liste_chr))+'\n')
		else:
			sys.exit('Orientation information is incorrect')
		l1 = file.readline().split()
		l2 = file.readline().split()
	outfile.close()
	file.close()


def trie_discord_rf(UN, DEUX, MINI_DIS, MINI, MAXI, LISTE_CHR):
	ident1 = UN[0]
	ident2 = DEUX[0]
	flag1 = UN[1]
	flag2 = DEUX[1]
	pos1 = UN[3]
	pos2 = DEUX[3]
	chr1 = UN[2]
	chr2 = DEUX[2]
	INSERT1 = abs(int(UN[8]))
	INSERT2 = abs(int(DEUX[8]))
	if INSERT1 != INSERT2:
		sys.exit('The two mate do not share same absolute insert size')
	if ident1 != ident2:
		sys.exit('The reads are not paired')
	elif (flag1 == '83' and flag2 == '163'):					#read rf avec R en mate 1
		if INSERT1 > MAXI:											##le read est discordant de deletion
			type = 'del'
		elif INSERT1 < MINI:										##le read est discordant d insertion
			type = 'ins'
		else:														##le read est concordant
			type = 'ok'
	elif (flag1 == '99' and flag2 == '147'): 					#read rf avec F en mate 1
		if INSERT1 > MAXI:											##le read est discordant de deletion
			type = 'del'
		elif INSERT1 < MINI:										##le read est discordant d insertion
			type = 'ins'
		else:														##le read est concordant
			type = 'ok'
	else:
		if chr1 != chr2:										#discordance de chromosome
			if (flag1 == '113' and flag2 == '177'):					##discordance "rr"
				type = 'chr_rr'
			elif (flag1 == '65' and flag2 == '129'):				##discordance "ff"
				type = 'chr_ff'
			elif (flag1 == '81' and flag2 == '161'):				##la mate 1 est Reverse la 2 est Forward
				if LISTE_CHR.index(chr1) < LISTE_CHR.index(chr2):
					type = 'chr_rf'										###Le chr mate1 avant chr mate2
				else:
					type = 'chr_fr'										###Le chr mate1 apres chr mate2
			elif (flag1 == '97' and flag2 == '145'):				##la mate 1 est Forward la 2 est Reverse
				if LISTE_CHR.index(chr1) < LISTE_CHR.index(chr2):
					type = 'chr_fr'										###Le chr mate1 avant chr mate2
				else:
					type = 'chr_rf'										###Le chr mate1 apres chr mate2
			else:
				sys.exit('Problem in sam flag')
		elif (flag1 == '113' and flag2 == '177'):				#discordance rr
			type = 'rr'
		elif (flag1 == '65' and flag2 == '129'):				#discordance ff
			type = 'ff'
		elif (flag1 == '81' and flag2 == '161'):				#read rf ou fr et la mate 1 est Reverse
			if int(pos1) <= int(pos2):								##read rf
				if INSERT1 > int(MAXI):									###le read est discordant de deletion
					type = 'del'
				elif INSERT1 < int(MINI):								###le read est discordant d insertion
					type = 'ins'
				else:													###le read est concordant
					type = 'ok'
			else:													##read fr
				type = 'fr'
		elif (flag1 == '97' and flag2 == '145'):				#read rf ou fr et la mate 1 est Forward
			if int(pos1) >= int(pos2):								##read rf
				if INSERT1 > int(MAXI):									###le read est discordant de deletion
					type = 'del'
				elif INSERT1 < int(MINI):								###le read est discordant d insertion
					type = 'ins'
				else:													###le read est concordant
					type = 'ok'
			else:													##read fr
				type = 'fr'
		else:
			sys.exit('Problem in sam flag')
	if type in ['ff', 'rr', 'fr', 'del']:
		if INSERT1 < MINI_DIS:
			return [ident1, pos1, pos2, chr1, chr2, 'discard', type]
		else:
			return [ident1, pos1, pos2, chr1, chr2, type]
	else:
		return [ident1, pos1, pos2, chr1, chr2, type]

def trie_discord_fr(UN, DEUX, MINI_DIS, MINI, MAXI, LISTE_CHR):
	ident1 = UN[0]
	ident2 = DEUX[0]
	flag1 = UN[1]
	flag2 = DEUX[1]
	pos1 = UN[3]
	pos2 = DEUX[3]
	chr1 = UN[2]
	chr2 = DEUX[2]
	INSERT1 = abs(int(UN[8]))
	INSERT2 = abs(int(DEUX[8]))
	if INSERT1 != INSERT2:
		sys.exit('The two mate do not share same absolute insert size')
	if ident1 != ident2:
		sys.exit('The reads are not paired')
	elif (flag1 == '83' and flag2 == '163'):					#read fr avec R en mate 1
		if INSERT1 > MAXI:											##le read est discordant de deletion
			type = 'del'
		elif INSERT1 < MINI:										##le read est discordant d insertion
			type = 'ins'
		else:														##le read est concordant
			type = 'ok'
	elif (flag1 == '99' and flag2 == '147'): 					#read fr avec F en mate 1
		if INSERT1 > MAXI:											##le read est discordant de deletion
			type = 'del'
		elif INSERT1 < MINI:										##le read est discordant d insertion
			type = 'ins'
		else:														##le read est concordant
			type = 'ok'
	else:
		if chr1 != chr2:										#discordance de chromosome
			if (flag1 == '113' and flag2 == '177'):					##discordance "rr"
				type = 'chr_rr'
			elif (flag1 == '65' and flag2 == '129'):				##discordance "ff"
				type = 'chr_ff'
			elif (flag1 == '81' and flag2 == '161'):				##la mate 1 est Reverse la 2 est Forward
				if LISTE_CHR.index(chr1) < LISTE_CHR.index(chr2):
					type = 'chr_rf'										###Le chr mate1 avant chr mate2
				else:
					type = 'chr_fr'										###Le chr mate1 apres chr mate2
			elif (flag1 == '97' and flag2 == '145'):				##la mate 1 est Forward la 2 est Reverse
				if LISTE_CHR.index(chr1) < LISTE_CHR.index(chr2):
					type = 'chr_fr'										###Le chr mate1 avant chr mate2
				else:
					type = 'chr_rf'										###Le chr mate1 apres chr mate2
			else:
				sys.exit('Problem in sam flag')
		elif (flag1 == '113' and flag2 == '177'):				#discordance rr
			type = 'rr'
		elif (flag1 == '65' and flag2 == '129'):				#discordance ff
			type = 'ff'
		elif (flag1 == '81' and flag2 == '161'):				#read rf ou fr et la mate 1 est Reverse
			if int(pos1) >= int(pos2):								##read rf
				if INSERT1 > int(MAXI):									###le read est discordant de deletion
					type = 'del'
				elif INSERT1 < int(MINI):								###le read est discordant d insertion
					type = 'ins'
				else:													###le read est concordant
					type = 'ok'
			else:													##read fr
				type = 'rf'
		elif (flag1 == '97' and flag2 == '145'):				#read rf ou fr et la mate 1 est Forward
			if int(pos1) <= int(pos2):								##read rf
				if INSERT1 > int(MAXI):									###le read est discordant de deletion
					type = 'del'
				elif INSERT1 < int(MINI):								###le read est discordant d insertion
					type = 'ins'
				else:													###le read est concordant
					type = 'ok'
			else:													##read fr
				type = 'rf'
		else:
			sys.exit('Problem in sam flag')
	if type in ['ff', 'rr', 'rf', 'del']:
		if INSERT1 < MINI_DIS:
			return [ident1, pos1, pos2, chr1, chr2, 'discard', type]
		else:
			return [ident1, pos1, pos2, chr1, chr2, type]
	else:
		return [ident1, pos1, pos2, chr1, chr2, type]

def grep(FILE, OUT, MOT, KEEP, COL):
	outfile = open(OUT, 'w')
	i = 0
	if COL == 'all':
		if KEEP == 'keep':
			for line in open(FILE):
				if MOT in line.split():
					outfile.write(line)
					i += 1
		elif KEEP == 'rm':
			for line in open(FILE):
				if not(MOT in line.split()):
					outfile.write(line)
					i += 1
		else:
			sys.exit('bug in grep')
	else:
		if KEEP == 'keep':
			for line in open(FILE):
				if MOT in line.split():
					data = line.split()
					liste = []
					for n in COL:
						liste.append(data[n])
					if MOT in liste:
						outfile.write(line)
						i += 1
		elif KEEP == 'rm':
			for line in open(FILE):
				if MOT in line.split():
					data = line.split()
					liste = []
					for n in COL:
						liste.append(data[n])
					if not(MOT in liste):
						outfile.write(line)
						i += 1
				else:
					outfile.write(line)
					i += 1
		else:
			sys.exit('bug in grep')
	outfile.close()
	return i

def calcul_prop(OUT, SAM, INSERT, INFO_CHR):
	############################################
	#recording chromosome order
	############################################
	DIC_REF = {}
	file = open(INFO_CHR)
	for line in file:
		data = line.split()
		if data:
			DIC_REF[data[0]] = int(data[1])
	file.close()

	############################################
	#Calculating
	############################################
	OUTFILE = open(OUT,'w')
	FEN = int(INSERT)
	#now for each window size
	FILE  = open(SAM)
	NB = 0
	WINDOW = 1
	CHR = ''
	LISTE_CHR = []
	for LINE in FILE:
		DATA = LINE.split()
		if DATA != []:
			if LINE[0] != '@': # pour gerer les headers
				if CHR == '':
					CHR = DATA[2]
					LISTE_CHR.append(DATA[2])
				elif DATA[2] != CHR:#changing of chromosomes
					if WINDOW+FEN < DIC_REF[CHR]:
						OUTFILE.write(CHR+'\t'+str(WINDOW)+'\t'+str(WINDOW+FEN-1)+'\t'+str(NB)+'\n')
						WINDOW = WINDOW + FEN
						while DIC_REF[CHR] > (WINDOW + FEN):
							OUTFILE.write(CHR+'\t'+str(WINDOW)+'\t'+str(WINDOW+FEN-1)+'\t0\n')
							WINDOW = WINDOW + FEN
						OUTFILE.write(CHR+'\t'+str(WINDOW)+'\t'+str(DIC_REF[CHR])+'\t0\n')
					elif WINDOW >= DIC_REF[CHR]:
						#The window is larger than the sequence
						OUTFILE.write(CHR+'\t0\t'+str(DIC_REF[CHR])+'\t0\n')
						# sys.exit('ERROR : Probleme in the script, this is a crash')
					else:
						OUTFILE.write(CHR+'\t'+str(WINDOW)+'\t'+str(DIC_REF[CHR])+'\t'+str(NB)+'\n')
					CHR = DATA[2]
					NB = 0
					WINDOW = 1
					LISTE_CHR.append(DATA[2])
				if DATA[2] == CHR:
					if int(DATA[3]) < (WINDOW + FEN):
						NB = NB + 1
					else:
						if NB:
							OUTFILE.write(CHR+'\t'+str(WINDOW)+'\t'+str(WINDOW+FEN-1)+'\t'+str(NB)+'\n')
						elif WINDOW == 1:
							OUTFILE.write(CHR+'\t'+str(WINDOW)+'\t'+str(WINDOW+FEN-1)+'\t0\n')
						NB = 1
						WINDOW = WINDOW + FEN
						while int(DATA[3]) >= (WINDOW + FEN):
							OUTFILE.write(CHR+'\t'+str(WINDOW)+'\t'+str(WINDOW+FEN-1)+'\t0\n')
							WINDOW = WINDOW + FEN
	if WINDOW+FEN < DIC_REF[CHR]:
		OUTFILE.write(CHR+'\t'+str(WINDOW)+'\t'+str(WINDOW+FEN-1)+'\t'+str(NB)+'\n')
		WINDOW = WINDOW + FEN
		while DIC_REF[CHR] > (WINDOW + FEN):
			OUTFILE.write(CHR+'\t'+str(WINDOW)+'\t'+str(WINDOW+FEN-1)+'\t0\n')
			WINDOW = WINDOW + FEN
		OUTFILE.write(CHR+'\t'+str(WINDOW)+'\t'+str(DIC_REF[CHR])+'\t0\n')
	elif WINDOW >= DIC_REF[CHR]:
		#The window is larger than the sequence
		OUTFILE.write(CHR+'\t0\t'+str(DIC_REF[CHR])+'\t0\n')
		# sys.exit('ERROR : Probleme in the script, this is a crash')
	else:
		OUTFILE.write(CHR+'\t'+str(WINDOW)+'\t'+str(DIC_REF[CHR])+'\t'+str(NB)+'\n')
	for n in DIC_REF:
		if not(n in LISTE_CHR):
			WINDOW = 1
			if WINDOW+FEN < DIC_REF[n]:
				OUTFILE.write(n+'\t'+str(WINDOW)+'\t'+str(WINDOW+FEN-1)+'\t0\n')
				WINDOW = WINDOW + FEN
				while DIC_REF[n] > (WINDOW + FEN):
					OUTFILE.write(n+'\t'+str(WINDOW)+'\t'+str(WINDOW+FEN-1)+'\t0\n')
					WINDOW = WINDOW + FEN
				OUTFILE.write(n+'\t'+str(WINDOW)+'\t'+str(DIC_REF[n])+'\t0\n')
			elif WINDOW >= DIC_REF[n]:
				#The window is larger than the sequence
				OUTFILE.write(n+'\t0\t'+str(DIC_REF[n])+'\t0\n')
				# sys.exit('ERROR : Probleme in the script, this is a crash')
			else:
				OUTFILE.write(n+'\t'+str(WINDOW)+'\t'+str(DIC_REF[n])+'\t0\n')
	OUTFILE.close()

def prop_dis(N_DIS, DIS, OUT, INFO_CHR):
	############################################
	#recording chromosome order
	############################################
	DIC_REF = {}
	file = open(INFO_CHR)
	for line in file:
		data = line.split()
		if data:
			DIC_REF[data[0]] = int(data[1])
	file.close()

	outfile = open(OUT,'w')
	CHR = ''
	for n in DIC_REF:
		FILE = open(DIS)
		DATA = FILE.readline().split()
		if DATA == []:
			break
		while DATA[0] != n:
			DATA = FILE.readline().split()
			if DATA == []:
				break
		FILE_N = open(N_DIS)
		DATA_N = FILE_N.readline().split()
		if DATA_N == []:
			break
		while DATA_N[0] != n:
			DATA_N = FILE_N.readline().split()
			if DATA_N == []:
				break
		# print n, DATA_N, 'toto', DATA
		while DATA_N[0] == n:
			if DATA_N[1] == DATA[1] and DATA_N[2] == DATA[2]:
				if int(DATA_N[3]) == int(DATA[3]) == 0:
					outfile.write(DATA[0]+'\t'+DATA[1]+'\t'+DATA[2]+'\t'+str(-0.2)+'\n')
				else:
					outfile.write(DATA[0]+'\t'+DATA[1]+'\t'+DATA[2]+'\t'+str(int(DATA[3])/float(int(DATA_N[3])+int(DATA[3])))+'\n')
			else:
				sys.exit('ERROR : windows do not match')
			DATA_N = FILE_N.readline().split()
			DATA = FILE.readline().split()
			if DATA_N == []:
				break
		FILE_N.close()
		FILE.close()
	outfile.close()


def calcul_discord_prop_and_parse(LOCA_PROGRAMS, CHR, SAM, LISTE_TYPE, OUT_INS, OUT_DEL, OUT_FR, OUT_RF, OUT_FF, OUT_RR, OUT_CHR_FR, OUT_CHR_RF, OUT_CHR_FF, OUT_CHR_RR, OUT_DISCARDED, DISCORD_PROP, OR, EXCLUDE):
	to_exclude = EXCLUDE.split("=")
	LISTE_FILTERED = LISTE_TYPE+'.filtered'
	file = open(LISTE_TYPE)
	outfile = open(LISTE_FILTERED,'w')
	for line in file:
		data = line.split()
		if data:
			absent = 1
			for n in to_exclude:
				if n in data:
					absent = 0
			if absent:
				outfile.write(line)
	outfile.close()

	############################################
	#Calculating discordant proportions
	############################################
	liste_discord = (tempfile.NamedTemporaryFile().name)+'.list'
	grep(LISTE_FILTERED, liste_discord, 'ok', 'rm', [5])
	read_discord = (tempfile.NamedTemporaryFile().name)+'.sam'
	read_non_discord = (tempfile.NamedTemporaryFile().name)+'.sam'
	filter1 = '%s -jar %s FilterSamReads INPUT=%s OUTPUT=%s SORT_ORDER=coordinate QUIET=true FILTER=includeReadList READ_LIST_FILE=%s WRITE_READS_FILES=false MAX_RECORDS_IN_RAM=5000000 VERBOSITY=WARNING VALIDATION_STRINGENCY=SILENT' % (LOCA_PROGRAMS.get('Programs','java'), LOCA_PROGRAMS.get('Programs','picard-tool'), SAM, read_discord, liste_discord)
	filter2 = '%s -jar %s FilterSamReads INPUT=%s OUTPUT=%s SORT_ORDER=coordinate QUIET=true FILTER=excludeReadList READ_LIST_FILE=%s WRITE_READS_FILES=false MAX_RECORDS_IN_RAM=5000000 VERBOSITY=WARNING VALIDATION_STRINGENCY=SILENT' % (LOCA_PROGRAMS.get('Programs','java'), LOCA_PROGRAMS.get('Programs','picard-tool'), SAM, read_non_discord, liste_discord)
	run_job (filter1, 'Error in FilterSamReads.jar:')
	run_job (filter2, 'Error in FilterSamReads.jar:')

	discord_count = (tempfile.NamedTemporaryFile().name)+'.count'
	non_discord_count = (tempfile.NamedTemporaryFile().name)+'.count'
	calcul_prop(discord_count, read_discord, 1000, CHR)
	calcul_prop(non_discord_count, read_non_discord, 1000, CHR)
	prop_dis(non_discord_count, discord_count, DISCORD_PROP, CHR)


	# print liste_discord
	# print discord_count
	# print non_discord_count
	os.remove(discord_count)
	os.remove(non_discord_count)
	os.remove(liste_discord)

	############################################
	#Parsing bam file
	############################################

	LISTE_EMPTY = []

	if OR == 'rf':
		sam2bam = '%s view -bS %s > %s' % (LOCA_PROGRAMS.get('Programs','samtools'), read_non_discord ,OUT_RF)
		liste_discord_orient = (tempfile.NamedTemporaryFile().name)+'orient.list'
		nb_dis = grep(LISTE_FILTERED, liste_discord_orient, 'fr', 'keep', [5])
		filter1 = '%s -jar %s FilterSamReads INPUT=%s OUTPUT=%s SORT_ORDER=coordinate QUIET=true FILTER=includeReadList READ_LIST_FILE=%s WRITE_READS_FILES=false MAX_RECORDS_IN_RAM=5000000 VERBOSITY=WARNING VALIDATION_STRINGENCY=SILENT' % (LOCA_PROGRAMS.get('Programs','java'), LOCA_PROGRAMS.get('Programs','picard-tool'), SAM, OUT_FR, liste_discord_orient)
	elif OR == 'fr':
		sam2bam = '%s view -bS %s > %s' % (LOCA_PROGRAMS.get('Programs','samtools'), read_non_discord ,OUT_FR)
		liste_discord_orient = (tempfile.NamedTemporaryFile().name)+'orient.list'
		nb_dis = grep(LISTE_FILTERED, liste_discord_orient, 'rf', 'keep', [5])
		filter1 = '%s -jar %s FilterSamReads INPUT=%s OUTPUT=%s SORT_ORDER=coordinate QUIET=true FILTER=includeReadList READ_LIST_FILE=%s WRITE_READS_FILES=false MAX_RECORDS_IN_RAM=5000000 VERBOSITY=WARNING VALIDATION_STRINGENCY=SILENT' % (LOCA_PROGRAMS.get('Programs','java'), LOCA_PROGRAMS.get('Programs','picard-tool'), SAM, OUT_RF, liste_discord_orient)
	else:
		mot = 'Unrecognized --orient option : '+OR
		sys.exit(mot)

	run_job (sam2bam, 'Error in samtools view:')

	if nb_dis > 0:
		run_job (filter1, 'Error in FilterSamReads.jar:')
	else:
		if OR == 'rf':
			outfile=open(OUT_FR,'w')
			outfile.close()
			LISTE_EMPTY.append(OUT_FR)
		elif OR == 'fr':
			outfile=open(OUT_RF,'w')
			outfile.close()
			LISTE_EMPTY.append(OUT_RF)

	liste_discord_ins = (tempfile.NamedTemporaryFile().name)+'ins.list'
	nb_ins = grep(LISTE_FILTERED, liste_discord_ins, 'ins', 'keep', [5])
	filter_ins = '%s -jar %s FilterSamReads INPUT=%s OUTPUT=%s SORT_ORDER=coordinate QUIET=true FILTER=includeReadList READ_LIST_FILE=%s WRITE_READS_FILES=false MAX_RECORDS_IN_RAM=5000000 VERBOSITY=WARNING VALIDATION_STRINGENCY=SILENT' % (LOCA_PROGRAMS.get('Programs','java'), LOCA_PROGRAMS.get('Programs','picard-tool'), SAM, OUT_INS, liste_discord_ins)

	liste_discord_del = (tempfile.NamedTemporaryFile().name)+'del.list'
	nb_del = grep(LISTE_FILTERED, liste_discord_del, 'del', 'keep', [5])
	filter_del = '%s -jar %s FilterSamReads INPUT=%s OUTPUT=%s SORT_ORDER=coordinate QUIET=true FILTER=includeReadList READ_LIST_FILE=%s WRITE_READS_FILES=false MAX_RECORDS_IN_RAM=5000000 VERBOSITY=WARNING VALIDATION_STRINGENCY=SILENT' % (LOCA_PROGRAMS.get('Programs','java'), LOCA_PROGRAMS.get('Programs','picard-tool'), SAM, OUT_DEL, liste_discord_del)

	liste_discord_ff = (tempfile.NamedTemporaryFile().name)+'ff.list'
	nb_ff = grep(LISTE_FILTERED, liste_discord_ff, 'ff', 'keep', [5])
	filter_ff = '%s -jar %s FilterSamReads INPUT=%s OUTPUT=%s SORT_ORDER=coordinate QUIET=true FILTER=includeReadList READ_LIST_FILE=%s WRITE_READS_FILES=false MAX_RECORDS_IN_RAM=5000000 VERBOSITY=WARNING VALIDATION_STRINGENCY=SILENT' % (LOCA_PROGRAMS.get('Programs','java'), LOCA_PROGRAMS.get('Programs','picard-tool'), SAM, OUT_FF, liste_discord_ff)

	liste_discord_rr = (tempfile.NamedTemporaryFile().name)+'rr.list'
	nb_rr = grep(LISTE_FILTERED, liste_discord_rr, 'rr', 'keep', [5])
	filter_rr = '%s -jar %s FilterSamReads INPUT=%s OUTPUT=%s SORT_ORDER=coordinate QUIET=true FILTER=includeReadList READ_LIST_FILE=%s WRITE_READS_FILES=false MAX_RECORDS_IN_RAM=5000000 VERBOSITY=WARNING VALIDATION_STRINGENCY=SILENT' % (LOCA_PROGRAMS.get('Programs','java'), LOCA_PROGRAMS.get('Programs','picard-tool'), SAM, OUT_RR, liste_discord_rr)

	liste_discord_chr_ff = (tempfile.NamedTemporaryFile().name)+'chr_ff.list'
	nb_chr_ff = grep(LISTE_FILTERED, liste_discord_chr_ff, 'chr_ff', 'keep', [5])
	filter_chr_ff = '%s -jar %s FilterSamReads INPUT=%s OUTPUT=%s SORT_ORDER=coordinate QUIET=true FILTER=includeReadList READ_LIST_FILE=%s WRITE_READS_FILES=false MAX_RECORDS_IN_RAM=5000000 VERBOSITY=WARNING VALIDATION_STRINGENCY=SILENT' % (LOCA_PROGRAMS.get('Programs','java'), LOCA_PROGRAMS.get('Programs','picard-tool'), SAM, OUT_CHR_FF, liste_discord_chr_ff)

	liste_discord_chr_rr = (tempfile.NamedTemporaryFile().name)+'chr_rr.list'
	nb_chr_rr = grep(LISTE_FILTERED, liste_discord_chr_rr, 'chr_rr', 'keep', [5])
	filter_chr_rr = '%s -jar %s FilterSamReads INPUT=%s OUTPUT=%s SORT_ORDER=coordinate QUIET=true FILTER=includeReadList READ_LIST_FILE=%s WRITE_READS_FILES=false MAX_RECORDS_IN_RAM=5000000 VERBOSITY=WARNING VALIDATION_STRINGENCY=SILENT' % (LOCA_PROGRAMS.get('Programs','java'), LOCA_PROGRAMS.get('Programs','picard-tool'), SAM, OUT_CHR_RR, liste_discord_chr_rr)

	liste_discord_chr_fr = (tempfile.NamedTemporaryFile().name)+'chr_fr.list'
	nb_chr_fr = grep(LISTE_FILTERED, liste_discord_chr_fr, 'chr_fr', 'keep', [5])
	filter_chr_fr = '%s -jar %s FilterSamReads INPUT=%s OUTPUT=%s SORT_ORDER=coordinate QUIET=true FILTER=includeReadList READ_LIST_FILE=%s WRITE_READS_FILES=false MAX_RECORDS_IN_RAM=5000000 VERBOSITY=WARNING VALIDATION_STRINGENCY=SILENT' % (LOCA_PROGRAMS.get('Programs','java'), LOCA_PROGRAMS.get('Programs','picard-tool'), SAM, OUT_CHR_FR, liste_discord_chr_fr)

	liste_discord_chr_rf = (tempfile.NamedTemporaryFile().name)+'chr_rf.list'
	nb_chr_rf = grep(LISTE_FILTERED, liste_discord_chr_rf, 'chr_rf', 'keep', [5])
	filter_chr_rf = '%s -jar %s FilterSamReads INPUT=%s OUTPUT=%s SORT_ORDER=coordinate QUIET=true FILTER=includeReadList READ_LIST_FILE=%s WRITE_READS_FILES=false MAX_RECORDS_IN_RAM=5000000 VERBOSITY=WARNING VALIDATION_STRINGENCY=SILENT' % (LOCA_PROGRAMS.get('Programs','java'), LOCA_PROGRAMS.get('Programs','picard-tool'), SAM, OUT_CHR_RF, liste_discord_chr_rf)

	liste_discord_discarded = (tempfile.NamedTemporaryFile().name)+'discarded.list'
	nb_discarded = grep(LISTE_FILTERED, liste_discord_discarded, 'discard', 'keep', [5])
	filter_discarded = '%s -jar %s FilterSamReads INPUT=%s OUTPUT=%s SORT_ORDER=coordinate QUIET=true FILTER=includeReadList READ_LIST_FILE=%s WRITE_READS_FILES=false MAX_RECORDS_IN_RAM=5000000 VERBOSITY=WARNING VALIDATION_STRINGENCY=SILENT' % (LOCA_PROGRAMS.get('Programs','java'), LOCA_PROGRAMS.get('Programs','picard-tool'), SAM, OUT_DISCARDED, liste_discord_discarded)

	if nb_ins > 0:
		run_job (filter_ins, 'Error in FilterSamReads.jar:')
	else:
		outfile=open(OUT_INS,'w')
		outfile.close()
		LISTE_EMPTY.append(OUT_INS)

	if nb_del > 0:
		run_job (filter_del, 'Error in FilterSamReads.jar:')
	else:
		outfile=open(OUT_DEL,'w')
		outfile.close()
		LISTE_EMPTY.append(OUT_DEL)

	if nb_ff > 0:
		run_job (filter_ff, 'Error in FilterSamReads.jar:')
	else:
		outfile=open(OUT_FF,'w')
		outfile.close()
		LISTE_EMPTY.append(OUT_FF)

	if nb_rr > 0:
		run_job (filter_rr, 'Error in FilterSamReads.jar:')
	else:
		outfile=open(OUT_RR,'w')
		outfile.close()
		LISTE_EMPTY.append(OUT_RR)

	if nb_chr_ff > 0:
		run_job (filter_chr_ff, 'Error in FilterSamReads.jar:')
	else:
		outfile=open(OUT_CHR_FF,'w')
		outfile.close()
		LISTE_EMPTY.append(OUT_CHR_FF)

	if nb_chr_rr > 0:
		run_job (filter_chr_rr, 'Error in FilterSamReads.jar:')
	else:
		outfile=open(OUT_CHR_RR,'w')
		outfile.close()
		LISTE_EMPTY.append(OUT_CHR_RR)

	if nb_chr_fr > 0:
		run_job (filter_chr_fr, 'Error in FilterSamReads.jar:')
	else:
		outfile=open(OUT_CHR_FR,'w')
		outfile.close()
		LISTE_EMPTY.append(OUT_CHR_FR)

	if nb_chr_rf > 0:
		run_job (filter_chr_rf, 'Error in FilterSamReads.jar:')
	else:
		outfile=open(OUT_CHR_RF,'w')
		outfile.close()
		LISTE_EMPTY.append(OUT_CHR_RF)

	if nb_discarded > 0:
		run_job (filter_discarded, 'Error in FilterSamReads.jar:')
	else:
		outfile=open(OUT_DISCARDED,'w')
		outfile.close()
		LISTE_EMPTY.append(OUT_DISCARDED)

	# print read_discord
	# print read_non_discord
	# print liste_discord_orient
	# print liste_discord_ins
	# print liste_discord_del
	# print liste_discord_ff
	# print liste_discord_rr
	# print liste_discord_chr_ff
	# print liste_discord_chr_rr
	# print liste_discord_chr_fr
	# print liste_discord_chr_rf
	os.remove(read_discord)
	os.remove(read_non_discord)
	os.remove(liste_discord_orient)
	os.remove(liste_discord_ins)
	os.remove(liste_discord_del)
	os.remove(liste_discord_ff)
	os.remove(liste_discord_rr)
	os.remove(liste_discord_chr_ff)
	os.remove(liste_discord_chr_rr)
	os.remove(liste_discord_chr_fr)
	os.remove(liste_discord_chr_rf)
	os.remove(liste_discord_discarded)
	return LISTE_EMPTY



def __main__():
	#Parse Command Line
	parser = optparse.OptionParser(usage="python %prog [options]\n\nProgram designed by Guillaume MARTIN : guillaume.martin@cirad.fr\n\nThis program take in input a sam/bam file,"
	" calculate proportion of discordant reads on 1kb window size and parse the sam/bam file in several bam files corresponding to the different discordant types")

	# Wrapper options.
	parser.add_option( '', '--sam', dest='sam', default='not_filled', help='Paired sam/bam file')
	parser.add_option( '', '--chr', dest='chr', default='not_filled', help='The tabulated file containing in col 1 : chromosome name, col 2: chromosome length. A line for each chromosomes')
	parser.add_option( '', '--sort', dest='sort', default='unsorted', help='The sort order of the sam/bam input: coordinate, queryname, unsorted')
	parser.add_option( '', '--type', dest='type', default='sam', help='Input type : sam or bam, [default: %default]')
	parser.add_option( '', '--rminput', dest='rminput', default='n', help='remove sam/bam input file: y or n, [default: %default]')
	parser.add_option( '', '--orient', dest='orient', default='rf', help='The expected orientation: rf or fr, [default: %default]')
	parser.add_option( '', '--mini_dis', dest='mini_dis', default='10000', help='The minimal insert size to keep the discordant read for structural variation search (integer), [default: %default]')
	parser.add_option( '', '--mini', dest='mini', default='not_filled', help='The minimum insert size (integer)')
	parser.add_option( '', '--maxi', dest='maxi', default='not_filled', help='The maximum insert size (integer)')
	parser.add_option( '', '--out_ins', dest='out_ins', default='discord_ins.bam', help='Output bam file of insertion discordant read, [default: %default]')
	parser.add_option( '', '--out_del', dest='out_del', default='discord_del.bam', help='Output bam file of deletion discordant read, [default: %default]')
	parser.add_option( '', '--out_fr', dest='out_fr', default='discord_fr.bam', help='Output bam file mapping in fr and not contained in "out_ins" of "out_del" if appropiate, [default: %default]')
	parser.add_option( '', '--out_rf', dest='out_rf', default='discord_rf.bam', help='Output bam file mapping in rf and not contained in "out_ins" of "out_del" if appropiate, [default: %default]')
	parser.add_option( '', '--out_ff', dest='out_ff', default='discord_ff.bam', help='Output bam file mapping in ff, [default: %default]')
	parser.add_option( '', '--out_rr', dest='out_rr', default='discord_rr.bam', help='Output bam file mapping in rr, [default: %default]')
	parser.add_option( '', '--out_chr_fr', dest='out_chr_fr', default='discord_chr_fr.bam', help='Output bam file mapping in chr_fr, [default: %default]')
	parser.add_option( '', '--out_chr_rf', dest='out_chr_rf', default='discord_chr_rf.bam', help='Output bam file mapping in chr_rf, [default: %default]')
	parser.add_option( '', '--out_chr_ff', dest='out_chr_ff', default='discord_chr_ff.bam', help='Output bam file mapping in chr_ff, [default: %default]')
	parser.add_option( '', '--out_chr_rr', dest='out_chr_rr', default='discord_chr_rr.bam', help='Output bam file mapping in chr_rr, [default: %default]')
	parser.add_option( '', '--out_discarded', dest='out_discarded', default='discarded.bam', help='Output bam file discarded by mini_dis parameter, [default: %default]')
	parser.add_option( '', '--liste_type', dest='liste_type', default='liste_type.txt', help='Output of liste of col 1: reads, col2: their discordant type, col3: position of mate1, col4: position of mate2, [default: %default]')
	parser.add_option( '', '--discord_prop', dest='discord_prop', default='discord_prop.txt', help='Output discordant read proportion on window of 1kb, [default: %default]')
	parser.add_option( '', '--exclude_chrom', dest='exclude_chrom', default='no_exclude', help='Exclude chromosomes from analysis. "no_exclude" or chromosomes names separated by "=", [default: %default]')
	parser.add_option( '', '--config', dest='config', default=None)
	(options, args) = parser.parse_args()



	pathname = os.path.dirname(sys.argv[0])

	loca_programs = ConfigParser.RawConfigParser()
	loca_programs.read(pathname+'/loca_programs.conf')

	if options.config:
		config = ConfigParser.RawConfigParser()
		config.read(options.config)
		if config.get('General','restimate') == 'n':
			trie2discord_pair(loca_programs, config.get('Remove_dup','out'), config.get('Remove_dup','type'), config.get('Remove_dup','sort'), float(config.get('General','mini_dis')), float(config.get('General','mini')), float(config.get('General','maxi')), config.get('General','orient'), config.get('General','chr'), options.liste_type)
		else:
			trie2discord_pair(loca_programs, config.get('Remove_dup','out'), config.get('Remove_dup','type'), config.get('Remove_dup','sort'), float(config.get('General','mini_dis')), float(config.get('Calc_coverage','mini')), float(config.get('Calc_coverage','maxi')), config.get('General','orient'), config.get('General','chr'), options.liste_type)
		empty = calcul_discord_prop_and_parse(loca_programs, config.get('General','chr'), config.get('Remove_dup','out'), options.liste_type, options.out_ins, options.out_del, options.out_fr, options.out_rf, options.out_ff, options.out_rr, options.out_chr_fr, options.out_chr_rf, options.out_chr_ff, options.out_chr_rr, options.out_discarded, options.discord_prop, config.get('General','orient'), config.get('General','exclude_chrom'))
		if config.get('Trie_discord','rminput') == 'y':
			os.remove(config.get('Remove_dup','out'))
		if options.out_ins in empty:
			config.set('Trie_discord', 'out_ins', 'empty')
		else:
			config.set('Trie_discord', 'out_ins', options.out_ins)
		if options.out_del in empty:
			config.set('Trie_discord', 'out_del', 'empty')
		else:
			config.set('Trie_discord', 'out_del', options.out_del)
		if options.out_fr in empty:
			config.set('Trie_discord', 'out_fr', 'empty')
		else:
			config.set('Trie_discord', 'out_fr', options.out_fr)
		if options.out_rf in empty:
			config.set('Trie_discord', 'out_rf', 'empty')
		else:
			config.set('Trie_discord', 'out_rf', options.out_rf)
		if options.out_ff in empty:
			config.set('Trie_discord', 'out_ff', 'empty')
		else:
			config.set('Trie_discord', 'out_ff', options.out_ff)
		if options.out_rr in empty:
			config.set('Trie_discord', 'out_rr', 'empty')
		else:
			config.set('Trie_discord', 'out_rr', options.out_rr)
		if options.out_chr_fr in empty:
			config.set('Trie_discord', 'out_chr_fr', 'empty')
		else:
			config.set('Trie_discord', 'out_chr_fr', options.out_chr_fr)
		if options.out_chr_rf in empty:
			config.set('Trie_discord', 'out_chr_rf', 'empty')
		else:
			config.set('Trie_discord', 'out_chr_rf', options.out_chr_rf)
		if options.out_chr_ff in empty:
			config.set('Trie_discord', 'out_chr_ff', 'empty')
		else:
			config.set('Trie_discord', 'out_chr_ff', options.out_chr_ff)
		if options.out_chr_rr in empty:
			config.set('Trie_discord', 'out_chr_rr', 'empty')
		else:
			config.set('Trie_discord', 'out_chr_rr', options.out_chr_rr)
		if options.out_discarded in empty:
			config.set('Trie_discord', 'out_discarded', 'empty')
		else:
			config.set('Trie_discord', 'out_discarded', options.out_discarded)
		config.set('Trie_discord', 'liste_type', options.liste_type)
		config.set('Trie_discord', 'discord_prop', options.discord_prop)
		config.set('Trie_discord', 'type', 'bam')
		with open(options.config, 'wb') as configfile:
			config.write(configfile)
	else:
		if options.sam == 'not_filled':
			mot = 'Please provide an argument for --sam'
			sys.exit(mot)
		if options.sort == 'not_filled':
			mot = 'Please provide an argument for --sort'
			sys.exit(mot)
		if options.chr == 'not_filled':
			mot = 'Please provide an argument for --chr'
			sys.exit(mot)
		if options.mini == 'not_filled':
			mot = 'Please provide an argument for --mini'
			sys.exit(mot)
		if options.maxi == 'not_filled':
			mot = 'Please provide an argument for --maxi'
			sys.exit(mot)
		trie2discord_pair(loca_programs, options.sam, options.type, options.sort, float(options.mini_dis), float(options.mini), float(options.maxi), options.orient, options.chr, options.liste_type)
		calcul_discord_prop_and_parse(loca_programs, options.chr, options.sam, options.liste_type, options.out_ins, options.out_del, options.out_fr, options.out_rf, options.out_ff, options.out_rr, options.out_chr_fr, options.out_chr_rf, options.out_chr_ff, options.out_chr_rr, options.out_discarded, options.discord_prop, options.orient, options.exclude_chrom)
		if options.rminput == 'y':
			os.remove(options.sam)
	os.remove(options.liste_type+".filtered")
if __name__ == "__main__": __main__()
