
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

import optparse, os, shutil, subprocess, sys, tempfile, fileinput, ConfigParser, operator, time, math, glob, datetime
from inspect import currentframe, getframeinfo

def stop_err( msg ):
	raise ValueError(msg)

def run_job (outLog, frameinfo, cmd_line, ERROR):
	logOutput = open(outLog, 'a')
	logOutput.write("\n"+cmd_line)
	logOutput.close()
	try:
		tmp = tempfile.NamedTemporaryFile().name
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
		stop_err( 'Line : '+str(frameinfo.lineno)+' - '+ERROR + str( e ) )

def run_job_silent (frameinfo, cmd_line, ERROR):
	try:
		tmp = tempfile.NamedTemporaryFile().name
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
		stop_err( 'Line : '+str(frameinfo.lineno)+'\n'+cmd_line+'\n'+ERROR + str( e ) )

def mediane(L):
	"""
		Give the median value of a list of integers

		:param L: A list of integers
		:type L: list
		:return: The median value
		:rtype: int
	"""

	L.sort()
	N = len(L)
	n = N/2.0
	p = int(n)
	if n == 0:
		return 0
	elif n == 1:
		return (L[0])
	elif n == p:
		return (L[p-1]+L[p])/2.0
	else:
		return float(L[p])

def extractSamFromPosition(LOCA_PROGRAMS, SAM, TYPE, CHR, START, END, OUT):
	"""
		Provide a sam file or bam file from coordinates

		This function create a sam or bam file, according to the input format, from a first sam or bam file and coordinates like chr01:10000-20000.

		:param LOCA_PROGRAMS: From the Configparser module. Contains the path of each programs
		:param SAM: The input sam or bam file
		:type SAM: str
		:param TYPE: The format of the input file. If the input file is a bam file, he must be indexed.
		:type TYPE: str ("sam" | "bam")
		:param CHR: File containing col1: chromosome name and col2: chromosome size
		:type CHR: str
		:param START: The start position
		:type START: int
		:param END: The end position
		:type END: int
		:param OUT: The name of the output file
		:return: void
	"""

	if TYPE == 'bam':
		bam2subbam = '%s view -bh %s %s:%s-%s -o %s' % (LOCA_PROGRAMS.get('Programs','samtools'), SAM, CHR, START, END, OUT)
		run_job_silent(getframeinfo(currentframe()), bam2subbam, 'Error in bam2subbam:\n')
	elif TYPE == 'sam':
		outbed = open('outbed_tmp', 'w')
		outbed.write(CHR+'\t'+START+'\t'+END)
		outbed.close()

		bam2subbam = '%s view -Sh -L %s -o %s %s' % (LOCA_PROGRAMS.get('Programs','samtools'), 'outbed_tmp', SAM, OUT)
		run_job_silent(getframeinfo(currentframe()), bam2subbam, 'Error in bam2subbam:\n')
		os.remove('outbed_tmp')
	else:
		raise ValueError('File type is invalid.')


####################################################################################################
#				Parse amont and aval locations
####################################################################################################

def amont_aval(LOCA_PROGRAMS, CHR, SAM, TYPE, OUT):
	"""
		Create two sam files with separate alignments according to the position of the reads.

		Parse a sam or bam file witch contains pair end alignments and separate the reads according to their position.
		When the reads are aligned on two different chromosomes, the order of the chromosomes is given by the parameter CHR.

		:param LOCA_PROGRAMS: From the Configparser module. Contains the path of each programs
		:param CHR: File containing col1: chromosome name and col2: chromosome size
		:type CHR: str
		:param SAM: The input sam or bam file.
		:type SAM: str
		:param TYPE: The format of the input file. If the input file is a bam file, he must be indexed.
		:type TYPE: str ("sam" | "bam")
		:param OUT: The name of the output file
		:type OUT: str
		:return: void

		.. seealso:: sort_sam()
		.. warnings:: This function only work in paired end reads.
	"""
	if TYPE == 'bam':
		bam2sam(LOCA_PROGRAMS, SAM, OUT+'_reverse2sam.sam')
		total = sort_sam(OUT+'_reverse2sam.sam', OUT, CHR)
		os.remove(OUT+'_reverse2sam.sam')
		return total
	elif TYPE == 'sam':
		total = sort_sam(SAM, OUT, CHR)
		return total
	else:
		raise ValueError('Unrecognised argument passed in --type option')


def sort_sam(SAM, OUT, CHR):

	"""
		Create two sam files with separate alignments according to the position of the reads.

		Parse a sam or bam file witch contains pair end alignments and separate the reads according to their position.
		When the reads are aligned on two different chromosomes, the order of the chromosomes is given by the parameter CHR.

		:param CHR: File containing col1: chromosome name and col2: chromosome size
		:type CHR: str
		:param SAM: The input sam or bam file.
		:type SAM: str
		:param OUT: The name of the output file
		:type OUT: str
		:return: The number of common reads in the two mate zones
		:rtype: int

		.. seealso:: amont_aval()
	"""

	#recording reference informations
	chrom_order = []

	listeZones = []
	file = open(CHR)
	for line in file:
		data = line.split()
		chrom_order.append(data[0])
	file.close()

	#parsing sam file
	amont = open(OUT+'_amont','w')
	aval = open(OUT+'_aval','w')
	file = open(SAM)
	total = 0
	for line in file:
		data = line.split()
		if data:
			if data[0][0] == '@': # lines corresponding to the sam header's
				amont.write(line)
				aval.write(line)
			else:
				if data[6] == '=' or data[6] == data[2]:
					if int(data[3]) <= int(data[7]):
						total += 1
						amont.write(line)
					elif int(data[3]) > int(data[7]):
						total += 1
						aval.write(line)
					else:
						raise ValueError('There seem to be a bug in sort_sam function')
				else:
					if chrom_order.index(data[2]) < chrom_order.index(data[6]):
						amont.write(line)
					else:
						aval.write(line)
					total += 1

	amont.close()
	aval.close()
	return total/2

def trie_sam(SAM, OUT, CHR):
	outfile1 = open(OUT+'_amont','w')
	outfile2 = open(OUT+'_aval','w')
	outfile3 = open(OUT+'_chr','w')
	file = open(SAM)
	total = 0
	# amont = 0
	# aval = 0
	for line in file:
		data = line.split()
		if data:
			if data[0][0] == '@':
				outfile1.write(line)
				outfile2.write(line)
				outfile3.write(line)
			else:
				if data[2] != CHR:
					total += 1
					outfile3.write(line)
				elif int(data[3]) < int(data[7]):
					total += 1
					outfile1.write(line)
				elif int(data[3]) > int(data[7]):
					total += 1
					outfile2.write(line)
	outfile1.close()
	outfile2.close()
	outfile3.close()
	return total/2


def triZonesAmontAval(zonesFile, CHR, OUT):
	"""
		Sort the two mate zones according to their position

		Sort the two mate zone identified by look_4_mate() according to the position of the first read.
		When the zones are on two different chromosomes, the order is defined if the CHR file.

		:param zonesFile: A file with 6 columns. First 3 columns corresponds to the coordinate of a zone (chr pos_start pos_end) and the 3 last columns corresponds to the coordinate of the mate zone.
		:type zonesFile: str
		:param CHR: File containing col1: chromosome name and col2: chromosome size
		:type CHR: str
		:param OUT: The name of the output file
		:type OUT: str
		:return: void

		.. seealso:: look_4_mate()
	"""

	chrom_order = []
	file = open(CHR)
	for line in file:
		words = line.split()
		if words:
			chrom_order.append(words[0])
	file.close()

	zones = open(zonesFile, 'r')
	output = open(OUT,'w')
	for line in zones:
		words = line.split()
		if words:
			if chrom_order.index(words[0]) < chrom_order.index(words[3]):
				output.write(line)
			elif chrom_order.index(words[0]) > chrom_order.index(words[3]):
				output.write('\t'.join(map(str, words[3:]))+'\t'+'\t'.join(map(str, words[0:3]))+'\n')
			elif chrom_order.index(words[0]) == chrom_order.index(words[3]):
				if int(words[1]) < int(words[4]):
					output.write(line)
				else:
					output.write('\t'.join(map(str, words[3:]))+'\t'+'\t'.join(map(str, words[0:3]))+'\n')
			else:
				raise ValueError('There is a problem in the chromosomes names')

####################################################################################################
#				For calculating discordant coverage
####################################################################################################
def calcul_cov(LOCA_PROGRAMS, SAM, TYPE, OUT):
	"""
		Calculate the coverage of a sam or bam file, site by site

		:param LOCA_PROGRAMS: From the Configparser module. Contains the path of each programs
		:param SAM: The input sam or bam file.
		:type SAM: str
		:param TYPE: The format of the input file
		:type TYPE: str ("sam" | "bam")
		:param OUT: The name of the output file
		:type OUT: str
		:return: void
	"""

	if TYPE == 'sam':
		#Convert the sam file to the bam format
		sam2bam = '%s view -bSh %s -o %s' % (LOCA_PROGRAMS.get('Programs','samtools'), SAM, SAM+'_BAM.bam')
		run_job_silent(getframeinfo(currentframe()), sam2bam, 'Error in sam2bam (calcul_cov):\n')
		SAM = SAM+'_BAM.bam'
	elif TYPE != 'bam':
		mot = TYPE+' argument passed in --type is not recognized'
		raise ValueError(mot)

	cal_cov = ('%s depth %s > %s') % (LOCA_PROGRAMS.get('Programs','samtools'), SAM, OUT)
	run_job_silent(getframeinfo(currentframe()), cal_cov, 'Error in calculating coverage:\n')

	#remove the intermediate bam file.
	if TYPE == 'sam':
		os.remove(SAM)


####################################################################################################
#				For zone identification
####################################################################################################


def select_sur_couv(COV, ZCOV, MAXCOV, MINCOV, NCOV, OUT):

	"""
		Detect zones according to the coverage and the contiguous non covered site

		:param COV: a tabular file with the coverage of each site
		:type COV: str
		:param ZCOV: The minimal number of covered site by zone
		:type ZCOV:
		:param MAXCOV: The maximum median coverage allowed
		:type MAXCOV: float
		:param MINCOV: The minimum median coverage allowed
		:type MINCOV: float
		:param NCOV: the number of contiguous uncovered sites allowed
		:type NCOV: int
		:param OUT: The name of the output file
		:type OUT: str
		:return: A list with the number of zones for the first element, and the length of the larger zone in the second element.
		:rtype: list
		.. seealso:: calcul_cov()
	"""

	debut, av = 0, 0
	chr = ""
	liste = []
	longerZone = 0
	outfile = open(OUT,'w')
	file = open(COV)
	i = 0
	for line in file:
		data = line.split()
		if data:
			POS = int(data[1])
			COV = int(data[2])
			if chr == '':  # on initialise sur le premier chromosome
				chr = data[0]
				debut = POS
				av = POS
				liste.append(COV)
			elif chr != data[0]:  # on change de chromosome
				val_med = mediane(liste)
				if len(liste) >= int(ZCOV) and val_med <= float(MAXCOV) and val_med >= float(MINCOV):
					outfile.write('\t'.join([chr,str(debut),str(av),str(len(liste)),str(val_med)])+'\n')

					if abs(int(av) - int(debut)) > longerZone:
						longerZone = abs(int(av) - int(debut))

					i += 1
				chr = data[0]
				debut = POS
				av = POS
				liste = []
				liste.append(COV)
			elif (POS - av) > int(NCOV):  # on change de zone, distance maximale pour une zone ou la couverture chute a 0
				val_med = mediane(liste)
				if len(liste) >= int(ZCOV) and val_med <= float(MAXCOV) and val_med >= float(MINCOV):
					outfile.write('\t'.join([chr,str(debut),str(av),str(len(liste)),str(val_med)])+'\n')

					if abs(int(av) - int(debut)) > longerZone:
						longerZone = abs(int(av) - int(debut))

					i += 1
				debut = POS
				av = POS
				liste = []
				liste.append(COV)
			else:
				liste.append(COV)
				av = POS
	outfile.close()
	return [i, longerZone]

def calculCovFromMateZones(LOCA_PROGRAMS, bamAmont, bamAval, TYPE, CHR_Amont, START_Amont, END_amont, CHR_Aval, START_Aval, END_Aval):
	"""
		Detect zones according to the coverage and the contiguous non covered site

		:param LOCA_PROGRAMS: From the Configparser module. Contains the path of each programs
		:param bamAmont: The bam file containing the upstream alignments
		:type bamAmont: str
		:param bamAval: The bam file containing the downstream alignments
		:type bamAval: str
		:param TYPE: The type of input file
		:type TYPE: str ("bam" | "sam")
		:param CHR_Amont: chromosome of the upstream zone
		:type CHR_Amont: float
		:param START_Amont: The start position of the upstream zone
		:type START_Amont: int
		:param END_amont: The end position of the upstream zone
		:type END_amont: int
		:param CHR_Aval: chromosome of the downstream zone
		:type CHR_Aval: float
		:param START_Aval: The start position of the downstream zone
		:type START_Aval: int
		:param END_Aval: The end position of the downstream zone
		:type END_Aval: int
		:return: A liste with 3 elements : the coverage of the upstream zone, the coverage of the downstream zone, and the number of commons reads from the two zones.
		:rtype: list
	"""

	#create sam and bam files for each zones
	extractSamFromPosition(LOCA_PROGRAMS, bamAmont, TYPE, CHR_Amont, START_Amont, END_amont, 'bam_amont_tmp')
	extractSamFromPosition(LOCA_PROGRAMS, bamAval, TYPE, CHR_Aval, START_Aval, END_Aval, 'bam_aval_tmp')

	#extract the reads in each zones
	extractReadsFromZones = '%s view %s | cut -f1' % (LOCA_PROGRAMS.get('Programs','samtools'), 'bam_amont_tmp')
	p = subprocess.Popen(extractReadsFromZones,stdout=subprocess.PIPE,stderr=subprocess.PIPE, shell=True)
	readsAmont, errors = p.communicate()

	extractReadsFromZones = '%s view %s | cut -f1' % (LOCA_PROGRAMS.get('Programs','samtools'), 'bam_aval_tmp')
	p = subprocess.Popen(extractReadsFromZones,stdout=subprocess.PIPE,stderr=subprocess.PIPE, shell=True)
	readsAval, errors = p.communicate()

	#find the common reads in the two mate zone
	readsAmont = filter(None, readsAmont.split('\n'))
	readsAval = filter(None, readsAval.split('\n'))
	commonReads = list(set(readsAmont) & set(readsAval))

	fcommonReads = open('commonReadsInMateZones_tmp', 'w')
	for item in commonReads:
		fcommonReads.write("%s\n" % item)
	fcommonReads.close()

	#create two sam files for each zones with only the common reads
	recupHeader = '%s view -H %s -o %s' % (LOCA_PROGRAMS.get('Programs','samtools'), 'bam_amont_tmp', 'sam_amont_selected_tmp')
	run_job_silent(getframeinfo(currentframe()), recupHeader, 'Error in recupHeader:\n')
	extractReads = '%s -R %s %s >> %s' % (LOCA_PROGRAMS.get('Programs','bamgrepreads'), 'commonReadsInMateZones_tmp', 'bam_amont_tmp', 'sam_amont_selected_tmp')
	run_job_silent(getframeinfo(currentframe()), extractReads, 'Error in extractReads:\n')

	recupHeader = '%s view -H %s -o %s' % (LOCA_PROGRAMS.get('Programs','samtools'), 'bam_aval_tmp', 'sam_aval_selected_tmp')
	run_job_silent(getframeinfo(currentframe()), recupHeader, 'Error in recupHeader:\n')
	extractReads = '%s -R %s %s >> %s' % (LOCA_PROGRAMS.get('Programs','bamgrepreads'), 'commonReadsInMateZones_tmp', 'bam_aval_tmp', 'sam_aval_selected_tmp')
	run_job_silent(getframeinfo(currentframe()), extractReads, 'Error in extractReads:\n')

	# calcul of the coverage of this two sam file
	calcul_cov(LOCA_PROGRAMS, 'sam_amont_selected_tmp', 'sam', 'cov_amont_tmp')

	zoneCov = open('cov_amont_tmp', 'r')
	listeCov = []

	for line in zoneCov:
		cols = line.split()
		listeCov.append(int(cols[2]))
	zoneCov.close()

	medianCovAmont = mediane(listeCov)

	calcul_cov(LOCA_PROGRAMS, 'sam_aval_selected_tmp', 'sam', 'cov_aval_tmp')

	zoneCov = open('cov_aval_tmp', 'r')
	listeCov = []

	for line in zoneCov:
		cols = line.split()
		listeCov.append(int(cols[2]))
	zoneCov.close()

	medianCovAval = mediane(listeCov)

	for filename in glob.glob('*_tmp'):
		os.remove(filename)

	return [medianCovAmont, medianCovAval, len(commonReads)]


####################################################################################################
#				For mate zone identification
####################################################################################################
def decoup_chr_bam(LOCA_PROGRAMS, BAM, CHR, OUT, TYPE, overlap):
	"""
		Split a bam file by chromosome and with 1Mb windows, covered by 500kb.

		chr01 : 		 1Mb				   1Mb
				<-------------------><-------------------->
				          <--------------------><-------------------->
								    1Mb				      1Mb

		:param LOCA_PROGRAMS: From the Configparser module. Contains the path of each programs
		:param BAM: The bam file to split
		:type BAM: str
		:param CHR: File containing col1: chromosome name and col2: chromsome size
		:type CHR: str
		:param TYPE: The type of input file
		:type TYPE: str ("bam" | "sam")
		:param OUT: The output file name
		:type OUT: str
		:param longerZone: The length of the longer zone identified
		:type longerZone: int
		:param overlap: The overlap between two bam file.
		:type overlap: int
		:return: void
	"""

	file = open(CHR)
	liste_to_remove = []

	for line in file:
		data = line.split()
		if data:
			debut = 1
			while debut <= int(data[1]):
				tempo2 = tempfile.NamedTemporaryFile()

				if (debut + 999999) < int(data[1]):
					fin = debut + 999999
				else:
					fin = int(data[1])

				if TYPE == 'bam':
					bam2subbam = '%s view %s %s:%s-%s | cut -f 1 > %s' %(LOCA_PROGRAMS.get('Programs','samtools'), BAM, data[0], debut, fin, tempo2.name)
				elif TYPE == 'sam':

					# create a bed file
					tempo1 = tempfile.NamedTemporaryFile()
					tempo1.write(data[0]+'\t'+debut+'\t'+fin+'\n')
					tempo1.flush()

					bam2subbam = '%s view -S -L %s %s | cut -f 1 > %s' % (LOCA_PROGRAMS.get('Programs','samtools'), tempo1.name, BAM, tempo2)

					# Convert the sam file to the bam format
					sam2bam = '%s view -bSh %s -o %s' % (LOCA_PROGRAMS.get('Programs','samtools'), BAM, BAM+'_BAM.bam')
					run_job_silent(getframeinfo(currentframe()), sam2bam, 'Error in sam2bam :\n')

					# Remove the SAM file converted
					os.remove(BAM)
					BAM = BAM+'_BAM.bam'
				else:
					mot = TYPE+' argument passed in --type is not recognized'
					raise ValueError(mot)
				run_job_silent(getframeinfo(currentframe()), bam2subbam, 'Error in decoup_chr_bam:\n')
				tempo2.flush()

				# Here we have read mapping in the zone but we need the second mate...
				if os.path.getsize(tempo2.name) != 0:
					recupHeader = '%s view -H %s -o %s' % (LOCA_PROGRAMS.get('Programs','samtools'), BAM, OUT+'_'+data[0]+'_'+str(debut)+'-'+str(fin)+'.sam')
					run_job_silent(getframeinfo(currentframe()), recupHeader, 'Error in recupHeader:\n')
					extractReads = '%s -R %s %s >> %s' % (LOCA_PROGRAMS.get('Programs','bamgrepreads'), tempo2.name, BAM, OUT+'_'+data[0]+'_'+str(debut)+'-'+str(fin)+'.sam')
					run_job_silent(getframeinfo(currentframe()), extractReads, 'Error in extractReads:\n')

					# Convert the sam file to the bam format for Samtools and bamgrepreads
					sam2bam = '%s view -bSh %s -o %s' % (LOCA_PROGRAMS.get('Programs','samtools'), OUT+'_'+data[0]+'_'+str(debut)+'-'+str(fin)+'.sam', OUT+'_'+data[0]+'_'+str(debut)+'-'+str(fin)+'.bam')
					run_job_silent(getframeinfo(currentframe()), sam2bam, 'Error in sam2bam :\n')

					# Remove the SAM file converted
					os.remove(OUT+'_'+data[0]+'_'+str(debut)+'-'+str(fin)+'.sam')

					# Indexation of the bam file
					index_bam_file(LOCA_PROGRAMS, OUT+'_'+data[0]+'_'+str(debut)+'-'+str(fin)+'.bam')

				liste_to_remove.append(BAM+'_'+data[0]+'.bam')
				tempo2.close()
				debut = debut + 1000000-overlap

	return liste_to_remove

def index_bam_file(LOCA_PROGRAMS, BAM):
	"""
		Index a bam file.

		The name of the index file is : input_name.bai


		:param LOCA_PROGRAMS: From the Configparser module. Contains the path of each programs
		:param BAM: The input bam file.
		:return: void
	"""

	indexBamFile = '%s index %s' % (LOCA_PROGRAMS.get('Programs','samtools'), BAM)
	run_job_silent(getframeinfo(currentframe()), indexBamFile, 'Error in index_bam_file:\n')

def sam2bam(LOCA_PROGRAMS, SAM, OUT):
	sam2bam = '%s view -bSh %s -o %s' % (LOCA_PROGRAMS.get('Programs','samtools'), SAM, OUT)
	run_job_silent(getframeinfo(currentframe()), sam2bam, 'Error in sam2bam:\n')

def bam2sam(LOCA_PROGRAMS, BAM, OUT):
	bam2sam = '%s view -h %s -o %s' % (LOCA_PROGRAMS.get('Programs','samtools'), BAM, OUT)
	run_job_silent(getframeinfo(currentframe()), bam2sam, 'Error in bam conversion to sam:\n')

def create_sub_sam(outLog, LOCA_PROGRAMS, BAM, CHR, START, END, OUT):

	tempo2 = tempfile.NamedTemporaryFile()
	bam2subbam = '%s view %s %s:%s-%s | cut -f 1 > %s' % (LOCA_PROGRAMS.get('Programs','samtools'), BAM, CHR, START, END, tempo2.name)
	run_job_silent(getframeinfo(currentframe()), bam2subbam, 'Error in create_sub_sam:\n')

	if os.path.getsize(tempo2.name) != 0:
		recupHeader = '%s view -H %s -o %s' % (LOCA_PROGRAMS.get('Programs','samtools'), BAM, OUT)
		run_job_silent(getframeinfo(currentframe()), recupHeader, 'Error in recupHeader:\n')
		extractReads = '%s -R %s %s >> %s' % (LOCA_PROGRAMS.get('Programs','bamgrepreads'), tempo2.name, BAM, OUT)
		run_job_silent(getframeinfo(currentframe()), extractReads, 'Error in extractReads:\n')
	else:
		logOutput = open(outLog, 'a')
		logOutput.write("\ntempo2.name : "+str(tempo2.name))
		logOutput.write("\nbam file : "+str(BAM))
		logOutput.write("\nCHR : "+str(CHR))
		logOutput.write("\nSTART : "+str(START))
		logOutput.write("\nEND : "+str(END))
		logOutput.write("\nwarning there is a bug in create_sub_sam")
		logOutput.close()

def search_dest(LOCA_PROGRAMS, TARGET, DEST, DEST_CHR, ECART, CHR, START, END, ZCOV, MAXCOV, MINCOV, OUT, TOTAL):
	"""
		Search destination zone from a first discordant zone

		From a discordant zone identified, this function search for mate zones according to the coverture

		:param LOCA_PROGRAMS: From the Configparser module. Contains the path of each programs
		:param TARGET: A sam file containing the alignments of a discordant zone
		:type TARGET: str
		:param DEST: A sam file containing the alignments of a mate discordant zone
		:type SAM: str
		:param ECART: The value that will be added or substracted to keep a read as the same destination.
		:type ECART: int
		:param CHR: File containing col1: chromosome name and col2: chromsome size
		:type CHR: str
		:param START: the start position of the discordant zone
		:param START: int
		:param END: the end position of the discordant zone
		:param END: int
		:param ZCOV: Minimal number of covered sites in the zone
		:param ZCOV: int
		:param MAXCOV: The maximal median coverage accepted
		:param MAXCOV: float
		:param MINCOV: The minimal median coverage accepted
		:param MINCOV: float
		:param OUT: The name of the output file
		:param OUT: int
		:param TOTAL: The number of common reads in
		:param TOTAL: int
		:return: void

		.. seealso:: amont_aval() look_4_mate() recal_border()
	"""

	outfile = open(OUT, 'a')
	#loading TARGET information
	file = open(TARGET)
	dico = set()

	for line in file:
		data = line.split()
		if data:
			if data[0][0] != '@':
				dico.add(line)
	file.close()

	while dico:
		size_liste = 0
		liste = []
		# Verification si le read est dans la zone de depart
		while size_liste == 0 and dico:
			line1 = dico.pop()
			decoupe = line1.split()
			if int(decoupe[3]) <= END and int(decoupe[3]) >= START and decoupe[2] == CHR:
				liste.append(line1)
				center = int(decoupe[7])
				L_dest = [center]
				chr_dest = decoupe[6]
				chr_target = decoupe[2]
				value = (((END - START) + 1) + ECART)
				deb_zone = (center - value)
				fin_zone = (center + value)
			size_liste = len(liste)
		size_liste = 0  # To initlsialise the loop

		while size_liste != len(liste):
			size_liste = len(liste)
			for n in dico:
				decoupe = n.split()
				if decoupe[6] == chr_dest and int(decoupe[3]) <= END and int(decoupe[3]) >= START and decoupe[2] == CHR:  # The read is in the target zone and the destination is on the same chromosome

					if int(decoupe[7]) >= deb_zone and int(decoupe[7]) <= fin_zone:  # The read has a similar destination
						L_dest.append(int(decoupe[7]))
						med_val = mediane(L_dest)
						value = (((END - START) + 1) + ECART)
						deb_zone = (med_val - value)
						fin_zone = (med_val + value)
						liste.append(n)
			for n in liste:
				if n in dico:
					dico.remove(n)

		if liste:
			tempo = tempfile.NamedTemporaryFile()

			for n in liste:
				data = n.split()
				tempo.write(data[0]+'\n')
			tempo.flush()

			info_target = recalc_border(LOCA_PROGRAMS, tempo.name, TARGET, ZCOV, MAXCOV, MINCOV)

			if chr_dest == '=':
				info_dest = recalc_border(LOCA_PROGRAMS, tempo.name, DEST, ZCOV, MAXCOV, MINCOV)
			else:
				info_dest = recalc_border(LOCA_PROGRAMS, tempo.name, DEST_CHR, ZCOV, MAXCOV, MINCOV)

			tempo.close()

			if info_target[4] == 'PASS' and info_dest[4] == 'PASS':
				outfile.write('\t'.join([info_target[0], info_target[1], info_target[2], info_dest[0], info_dest[1], info_dest[2]])+'\n')

	outfile.close()

def recalc_border(LOCA_PROGRAMS, LISTE, SAM, ZCOV, MAXCOV, MINCOV):
	"""
		Calculate the borders of a zone, based on the coverture

		:param LOCA_PROGRAMS: From the Configparser module. Contains the path of each programs
		:param LISTE: A one column file with le list of the common reads of a zone
		:type LISTE: str
		:param SAM: A sam file containing the alignments of a discordant zone
		:type SAM: str
		:param ZCOV: Minimal number of covered sites in the zone
		:param ZCOV: int
		:param MAXCOV: The maximal median coverage accepted
		:param MAXCOV: float
		:param MINCOV: The minimal median coverage accepted
		:param MINCOV: float

		:return: a list with the coordinate of the zone, the median coverture and the information if the zone is valid or not.
		:rtype: list
		.. seealso:: search_dest() look_4_mate()
	"""

	#Convert the sam file to the bam format for Samtools and bamgrepreads
	sam2bam(LOCA_PROGRAMS, SAM, SAM+'_bam')

	#extract the header of the bam file
	recupHeader = '%s view -H %s -o %s' % (LOCA_PROGRAMS.get('Programs','samtools'), SAM+'_bam', SAM+'_sub.sam')
	run_job_silent(getframeinfo(currentframe()), recupHeader, 'Error in (recalc_border) recupHeader :\n')

	#concatenate the alignments specified to the header
	extractReads = '%s -R %s %s >> %s' % (LOCA_PROGRAMS.get('Programs','bamgrepreads'), LISTE, SAM+'_bam', SAM+'_sub.sam')
	run_job_silent(getframeinfo(currentframe()), extractReads, 'Error in extractReads:\n')

	#Verify if the sam file contains alignments
	grepNbAlignments = 'grep -c -v -m 1 "^@" '+SAM+'_sub.sam'
	p = subprocess.Popen(grepNbAlignments,stdout=subprocess.PIPE,stderr=subprocess.PIPE, shell=True)
	nbAlignments, errors = p.communicate()

	if int(nbAlignments) > 0:
		#convert the sam file to bam file
		sam2bam(LOCA_PROGRAMS, SAM+'_sub.sam', SAM+'_sub.bam')
		calcul_cov(LOCA_PROGRAMS, SAM+'_sub.bam', 'bam', SAM+'_sub_cov')

		chr = ""
		file = open(SAM+'_sub_cov')
		liste_cov = []

		for line in file:
			data = line.split()
			if data != []:
				liste_cov.append(int(data[2]))
				if chr == "":
					chr = data[0]
					min = data[1]
					max = data[1]
				elif chr == data[0]:
					max = data[1]
				else:
					raise ValueError('There is a bug in recalc_border : several chromosomes are found in the cov file')
		MEDIAN = mediane(liste_cov)
		if len(liste_cov) >= ZCOV and MEDIAN >= MINCOV and MEDIAN <= MAXCOV:
			return [chr, min, max, MEDIAN, 'PASS']
		else:
			return [chr, min, max, MEDIAN, 'NO_PASS']
	else:
		return [0,0,0,0,'NO_PASS']

def look_4_mate(outLog, LOCA_PROGRAMS, TYPE, BAM, CHR, OUT, ZONE, ECART, ZCOV, MAXCOV, MINCOV, MINGAP, NB_ZONE, longerZone, medianInsert):
	"""
		From a discordant zone identified, search for a mate zone

		:param LOCA_PROGRAMS: From the Configparser module. Contains the path of each programs
		:param TYPE: the type of the input file
		:type TYPE: str ("sam" | "bam")
		:param BAM: A sam or bam file containing the alignments of a discordant zone
		:type BAM: str
		:param CHR: File containing col1: chromosome name and col2: chromsome size
		:type CHR: str
		:param OUT: the output file name
		:type OUT: str
		:param ZONE: A tabular file containing the discordant zone identified (chrom_name start_pos end_pos zone_lenght median_coverture)
		:type ZONE: str
		:param ZCOV: Minimal number of covered sites in the zone
		:type ZCOV: int
		:param MAXCOV: The maximal median coverage accepted
		:type MAXCOV: float
		:param MINCOV: The minimal median coverage accepted
		:type MINCOV: float
		:param MINGAP: The maximal distance of contiguous uncovered sites
		:type MINGAP: int
		:param NB_ZONE: The total number of discordant zones identified
		:type NB_ZONE: int
		:param longerZone: The length of the longer zone
		:type longerZone: int
		:param medianInsert: The median length of the insert
		:type medianInsert: int
		:return: void

		.. seealso:: calcul_cov(), select_sur_couv(), search_dest(), recal_border(), calcul_cov(), select_sur_couv()
	"""
	#loading chromosome informations in a dictionnary
	dico_chr_info = {}
	file = open(CHR)
	for line in file:
		data = line.split()
		if data:
			dico_chr_info[data[0]] = int(data[1])
	file.close()
	#for speed increase creation of a sub_bam file for each chromosome
	tmp = tempfile.NamedTemporaryFile().name
	current_dir_path = os.getcwd()
	out_tmp = current_dir_path+'/'+tmp.split('/')[-1]

	#calculate the overlap between two bam file
	overlap = int(max(longerZone*2, medianInsert*3))

	t0 = datetime.datetime.now()
	decoup_chr_bam(LOCA_PROGRAMS, BAM, CHR, out_tmp, TYPE, overlap)
	logOutput = open(outLog, 'a')
	logOutput.write("\nsplit the bam file : "+str(datetime.datetime.now()-t0))
	logOutput.close()
	#Now it's time to work on each zone
	outfile = open(OUT,'w')
	outfile.close()
	file = open(ZONE)
	i = 0
	t_look4mate = datetime.datetime.now()
	####This part to manage the sub_bam files
	chromosome = ""
	for line in file:
		data = line.split()
		if data != []:
			pos_zone_debut = int(data[1])
			pos_zone_fin = int(data[2])
			if chromosome != data[0]:
				chromosome = data[0]
				pos_debut = 1
				pos_fin = pos_debut + 999999
				if pos_fin > dico_chr_info[data[0]]:
					pos_fin = dico_chr_info[data[0]]

				while not((pos_debut <= pos_zone_debut) and (pos_zone_fin <= pos_fin)):
					pos_debut = pos_debut + 1000000-overlap
					pos_fin = pos_debut + 999999
					if pos_fin > dico_chr_info[data[0]]:
						pos_fin = dico_chr_info[data[0]]
			else:
				while not((pos_debut <= pos_zone_debut) and (pos_zone_fin <= pos_fin)):
					pos_debut = pos_debut + 1000000-overlap
					pos_fin = pos_debut + 999999
					if pos_fin > dico_chr_info[data[0]]:
						pos_fin = dico_chr_info[data[0]]

			#creation of a sub bam containing read mapping in each putative discordant zone and their mate
			create_sub_sam(outLog, LOCA_PROGRAMS, out_tmp+'_'+data[0]+'_'+str(pos_debut)+'-'+str(pos_fin)+'.bam', data[0], data[1], data[2], out_tmp+'_tempo_liste_'+data[0]+'.sam')

			total = trie_sam(out_tmp+'_tempo_liste_'+data[0]+'.sam', out_tmp+'_sam_'+data[0], data[0])

			#search for destination of second mate and perform selection
			search_dest(LOCA_PROGRAMS, out_tmp+'_sam_'+data[0]+'_amont', out_tmp+'_sam_'+data[0]+'_aval', out_tmp+'_sam_'+data[0]+'_chr', float(ECART), data[0], int(data[1]), int(data[2]), ZCOV, MAXCOV, MINCOV, OUT, total)
			search_dest(LOCA_PROGRAMS, out_tmp+'_sam_'+data[0]+'_aval', out_tmp+'_sam_'+data[0]+'_amont', out_tmp+'_sam_'+data[0]+'_chr', float(ECART), data[0], int(data[1]), int(data[2]), ZCOV, MAXCOV, MINCOV, OUT, total)

			logOutput = open(outLog, 'a')
			i += 1
			if i % 100 == 0:
				logOutput.write("\nEstimation for 7_step_remaining time for"+str( BAM)+ " : "+str( i) +"done in"+str( datetime.datetime.now() - t_look4mate)+ ","+str(  NB_ZONE-i)+ "remaining")
			logOutput.close()
	for filename in glob.glob(out_tmp+'_*'):
		os.remove(filename)

####################################################################################################
#				Merging identical and similar zones
####################################################################################################

def regroupe_filtre(DATA, DATA_prec, MAX):
	"""
		Search a potential overlap between a couple of two mate zones and another couple of mate zones

		:param DATA: list of a first couple mate zones (chr1 start_pos1 end_pos1 chr2 start_pos2 end_pos2)
		:type DATA: list
		:param DATA_prec:list of a second couple of mate zones (chr1 start_pos1 end_pos1 chr2 start_pos2 end_pos2)
		:type DATA_prec: list
		:param MAX: Maximal distance between two discordant zone to merge
		:type MAX: int
		:return: a list with the second couple of mate zone,  and the information if this couple is overlapping the first couple of mate zones.
		:rtype: list
		.. seealso:: merge_zone()
	"""

	debut1 = str(min(int(DATA_prec[1]),int(DATA[1])))
	fin1 = str(max(int(DATA_prec[2]),int(DATA[2])))
	debut2 = str(min(int(DATA_prec[4]),int(DATA[4])))
	fin2 = str(max(int(DATA_prec[5]),int(DATA[5])))
	size1 = str((int(fin1)-int(debut1)) + 1)
	size2 = str((int(fin2)-int(debut2)) + 1)

	if (DATA[0] == DATA_prec[0]) and (DATA[3] == DATA_prec[3]):
		if int(DATA_prec[1]) <= int(DATA[1]):
			if ((int(DATA[1]) - int(DATA_prec[2])) < MAX) and ((int(DATA[1]) < int(DATA[4]) and int(DATA_prec[1]) < int(DATA_prec[4])) or (int(DATA[1]) > int(DATA[4]) and int(DATA_prec[1]) > int(DATA_prec[4]))):
				if int(DATA[5]) <= int(DATA_prec[4]):#la zone arrive avant la zone precedente : ssssss pppppp
					if ((int(DATA_prec[4])-int(DATA[5])) < MAX):#on est sur deux zones contigues non chevauchantes
						return [[DATA[0], debut1, fin1, DATA_prec[3], debut2, fin2, DATA[6]], 'found']
					else:#les zones ne sont pas contigues
						return [DATA_prec, 'not_found']
				elif int(DATA_prec[5]) <= int(DATA[4]):#la zone arrive apres la zone precedente : pppppp ssssss
					if ((int(DATA[4])-int(DATA_prec[5])) < MAX):#on est sur deux zones contigues non chevauchantes
						return [[DATA[0], debut1, fin1, DATA_prec[3], debut2, fin2, DATA[6]], 'found']
					else:#les zones ne sont pas contigues
						return [DATA_prec, 'not_found']
				elif int(DATA[5]) >= int(DATA_prec[4]) and int(DATA[5]) <= int(DATA_prec[5]) and int(DATA[4]) <= int(DATA_prec[4]): #zone precedente chevauche sur la suivante en mappant apres : sssspspspppppp
					return [[DATA[0], debut1, fin1, DATA_prec[3], debut2, fin2, DATA[6]], 'found']
				elif int(DATA[5]) >= int(DATA_prec[5]) and int(DATA[4]) <= int(DATA_prec[4]): #la zone precedente est inclue dans la suivante  sssspspspspssss
					return [[DATA[0], debut1, fin1, DATA_prec[3], debut2, fin2, DATA[6]], 'found']
				elif int(DATA_prec[5]) >= int(DATA[4]) and int(DATA_prec[5]) <= int(DATA[5]) and int(DATA[4]) >= int(DATA_prec[4]): #zone precedente chevauche sur la suivante en mappant avant : pppppspspspssssss
					return [[DATA[0], debut1, fin1, DATA_prec[3], debut2, fin2, DATA[6]], 'found']
				elif int(DATA_prec[5]) >= int(DATA[5]) and int(DATA_prec[4]) <= int(DATA[4]): #la zone suivante est inclue dans la suivante  ppppppspspspspppppp
					return [[DATA[0], debut1, fin1, DATA_prec[3], debut2, fin2, DATA[6]], 'found']
				else:
					return [DATA_prec, 'not_found']
			else:
				return [DATA_prec, 'not_found']
		else: #DATA mapp before DATA_prec
			if ((int(DATA_prec[1]) - int(DATA[2])) < MAX) and ((int(DATA[1]) < int(DATA[4]) and int(DATA_prec[1]) < int(DATA_prec[4])) or (int(DATA[1]) > int(DATA[4]) and int(DATA_prec[1]) > int(DATA_prec[4]))):
				if int(DATA[5]) <= int(DATA_prec[4]):#la zone arrive avant la zone precedente : ssssss pppppp
					if ((int(DATA_prec[4])-int(DATA[5])) < MAX):#on est sur deux zones contigues non chevauchantes
						return [[DATA[0], debut1, fin1, DATA_prec[3], debut2, fin2, DATA[6]], 'found']
					else:#les zones ne sont pas contigues
						return [DATA_prec, 'not_found']
				elif int(DATA_prec[5]) <= int(DATA[4]):#la zone arrive apres la zone precedente : pppppp ssssss
					if ((int(DATA[4])-int(DATA_prec[5])) < MAX):#on est sur deux zones contigues non chevauchantes
						return [[DATA[0], debut1, fin1, DATA_prec[3], debut2, fin2, DATA[6]], 'found']
					else:#les zones ne sont pas contigues
						return [DATA_prec, 'not_found']
				elif int(DATA[5]) >= int(DATA_prec[4]) and int(DATA[5]) <= int(DATA_prec[5]) and int(DATA[4]) <= int(DATA_prec[4]): #zone precedente chevauche sur la suivante en mappant apres : sssspspspppppp
					return [[DATA[0], debut1, fin1, DATA_prec[3], debut2, fin2, DATA[6]], 'found']
				elif int(DATA[5]) >= int(DATA_prec[5]) and int(DATA[4]) <= int(DATA_prec[4]): #la zone precedente est inclue dans la suivante  sssspspspspssss
					return [[DATA[0], debut1, fin1, DATA_prec[3], debut2, fin2, DATA[6]], 'found']
				elif int(DATA_prec[5]) >= int(DATA[4]) and int(DATA_prec[5]) <= int(DATA[5]) and int(DATA[4]) >= int(DATA_prec[4]): #zone precedente chevauche sur la suivante en mappant avant : pppppspspspssssss
					return [[DATA[0], debut1, fin1, DATA_prec[3], debut2, fin2, DATA[6]], 'found']
				elif int(DATA_prec[5]) >= int(DATA[5]) and int(DATA_prec[4]) <= int(DATA[4]): #la zone suivante est inclue dans la suivante  ppppppspspspspppppp
					return [[DATA[0], debut1, fin1, DATA_prec[3], debut2, fin2, DATA[6]], 'found']
				else:
					return [DATA_prec, 'not_found']
			else:
				return [DATA_prec, 'not_found']
	else:
		return [DATA_prec, 'not_found']

def merge_zone(LOCA_PROGRAMS, CHR, FILE, MAX, bamAmont, bamAval, TYPE, OUT):

	"""
		Search if a couple of two mate zones can be merged with another couple of a two mate zones.

		:param LOCA_PROGRAMS: From the Configparser module. Contains the path of each programs
		:param CHR: File containing col1: chromosome name and col2: chromsome size
		:type CHR: str
		:param FILE:File containing all of the couples of mate zones (chr1 start_pos1 end_pos1 chr2 start_pos2 end_pos2)
		:type FILE: str
		:param MAX: Maximal distance between two discordant zone to merge
		:type MAX: int
		:param bamAmont: A sam or bam file containing the alignments of the upstream zone
		:type bamAmont: str
		:param bamAval: A sam or bam file containing the alignments of the downstream zone
		:type bamAval: str
		:param OUT: The output file name
		:rtype: str
		:return: void

		.. seealso:: regroupe_filtre()
	"""

	#Sort the zone file like : first 5 cols corresponds to the uphill zone, and the 5 cols following corresponds to the downhill zone
	triZonesAmontAval(FILE, CHR, FILE+'_sorted')

	#recording each zone in a dictionary
	data_prec = ''
	file = open(FILE+'_sorted')
	outfile = open(OUT,'w')
	dic = set()
	i = 0
	len_dic = 0

	for line in file:
		data = line.split()
		if data:
			i = i + 1
			if not(i in dic):
				data_prec = line.split()
				data_prec.append(i)
				dic.add(i)
				while len(dic) != len_dic:
					len_dic = len(dic)
					file2 = open(FILE+'_sorted')
					j = 0
					for line2 in file2:
						data_check = line2.split()
						if data_check:
							j = j + 1
							if not(j in dic):
								data2 = line2.split()
								data2.append(j)

								found = regroupe_filtre(data2, data_prec, MAX)
								data_prec = found[0]

								#We try the other order of zone 'amont and aval exchange places'
								if found[1] == 'not_found':
									new_data2 = [data2[3],data2[4],data2[5],data2[0],data2[1],data2[2],data2[6]]

									found = regroupe_filtre(new_data2, data_prec, MAX)
									data_prec = found[0]
								elif found[1] != 'found':
									raise ValueError('bug')
								dic.add(data_prec[6])

				covs = calculCovFromMateZones(LOCA_PROGRAMS, bamAmont, bamAval, TYPE, data_prec[0], data_prec[1], data_prec[2], data_prec[3], data_prec[4], data_prec[5])

				outfile.write(str(data_prec[0])+' '+str(data_prec[1])+' '+str(data_prec[2])+' '+str(int(data_prec[2])-int(data_prec[1]))+' '+str(covs[0])+' '+str(data_prec[3])+' '+str(data_prec[4])+' '+str(data_prec[5])+' '+str(int(data_prec[5])-int(data_prec[4]))+' '+str(covs[1])+' - '+str(covs[2])+'\n')
	outfile.close()

#########################################################################################################################################################
#				Calculating score
#########################################################################################################################################################
def calculate_score(FILE, YIS, MIS, YIC, MIC, MIN, CHR, OUT):

	"""
		Calculate a score of mate zones and create the final output file.

		:param FILE: The file containing all of the mate zones merged
		:type FILE: str
		:param YIS: The Y-intercept of the linear function for zone size that will give the first component of product giving the score
		:type YIS: int
		:param MIS: The minimal zone size for which the first component of product giving the score will be maximal
		:type MIS: int
		:param YIC: The Y-intercept of the linear function for coverage that will give the second component of product giving the score
		:type YIC: int
		:param MIC: The minimal zone coverage for which the second component of product giving the score will be maximal
		:type MIC: int
		:param MIN: The minimal score for a discordant zone to be identified as passed
		:type MIN: int
		:param CHR: File containing col1: chromosome name and col2: chromosome size
		:type CHR: str
		:param OUT: The output file name
		:type OUT: str
		:return: void

	"""
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

	#parametres affine sur multiplicateur de taille
	b = float(YIS)
	x = float(MIS)
	a = (1-b)/x

	#parametres affine sur multiplicateur de couverture
	d = float(YIC)
	y = float(MIC)
	c = (1-d)/y

	outfile = open(OUT, 'w')
	file = open(FILE)
	for line in file:
		data = line.split()
		if data:
			#score de depart
			score = 100
			#multiplicateur du score sur la taille
			if ((int(data[3])+int(data[8]))/2.0) >= x:
				score = score*1
			else:
				score = score*(a*((int(data[3])+int(data[8]))/2.0)+b)
			#multiplicateur du score sur la couverture
			if ((float(data[4])+float(data[9]))/2.0) >= y:
				score = score*1
			else:
				score = score*(c*((float(data[4])+float(data[9]))/2.0)+d)
			###############
			#Ordering data#
			###############
			if data[0] == data[5]:
				if int(data[1]) <= int(data[6]):
					data_out = list(data)
				else:
					data_out = [data[5], data[6], data[7], data[8], data[9], data[0], data[1], data[2], data[3], data[4], data[10], data[11]]
			else:
				if liste_chr.index(data[0]) < liste_chr.index(data[5]):
					data_out = list(data)
				else:
					data_out = [data[5], data[6], data[7], data[8], data[9], data[0], data[1], data[2], data[3], data[4], data[10], data[11]]
			data_out.append(str(score))
			if score >= MIN:
				data_out.append('PASSED')
			else:
				data_out.append('NOT_PASSED')
			outfile.write('\t'.join(data_out)+'\n')
	outfile.close()


def __main__():

	t_start = datetime.datetime.now()
	#Parse Command Line
	parser = optparse.OptionParser(usage="python %prog [options]\n\nProgram designed by Guillaume MARTIN : guillaume.martin@cirad.fr\n\nThis program take in input a sam/bam file,"
	" and calculate coverage for covered sites. The bam should be coordinate sorted")
	# Wrapper options.
	parser.add_option( '', '--ref', dest='ref', default='not_filled', help='The multifasta reference file')
	parser.add_option( '', '--sam', dest='sam', default='not_filled', help='Paired sam/bam file')
	parser.add_option( '', '--type', dest='type', default='bam', help='Input type : sam or bam, [default: %default]')
	parser.add_option( '', '--min_zone', dest='min_zone', default=500, help='Minimal number of covered sites in the zone, [default: %default]')
	parser.add_option( '', '--maxcov', dest='maxcov', default=300, help='The maximal median coverage accepted (float), [default: %default]')
	parser.add_option( '', '--mincov', dest='mincov', default=0, help='The minimal median coverage accepted (float), [default: %default]')
	parser.add_option( '', '--min_gap', dest='min_gap', default=300, help='The maximal distance of contiguous uncovered sites, [default: %default]')
	parser.add_option( '', '--ecart', dest='ecart', default=2000, help='The value that will be added or substracted to keep a read as the same destination. Recommended 3*sd(insert)')
	parser.add_option( '', '--chr', dest='chr', default='not_filled', help='File containing col1: chromosome name and col2: chromsome size')
	parser.add_option( '', '--med_insert', dest='med_insert', default='not_filled', help='The median insert size')
	parser.add_option( '', '--median_coverage', dest='median_coverage', default='not_filled', help='The median coverage')
	parser.add_option( '', '--mult_max_cov', dest='mult_max_cov', default='10', help='multiplicator of median coverage for maximal median coverage to keep a zone (float), [default: %default]')
	parser.add_option( '', '--mult_min_cov', dest='mult_min_cov', default='0.25', help='multiplicator of median coverage for minimal median coverage to keep a zone (float),  [default: %default]')
	parser.add_option( '', '--max_dist_merge', dest='max_dist_merge', default=1000, help='Maximal distance between two discordant zone to merge, [default: %default]')
	parser.add_option( '', '--YiS', dest='YiS', default=0, help='The Y-intercept of the linear function for zone size that will give the first component of product giving the score (integer), [default: %default]')
	parser.add_option( '', '--MiS', dest='MiS', default=1000, help='The minimal zone size for which the first component of product giving the score will be maximal (integer), [default: %default]')
	parser.add_option( '', '--YiC', dest='YiC', default=0, help='The Y-intercept of the linear function for coverage that will give the second component of product giving the score (integer), [default: %default]')
	parser.add_option( '', '--MiC', dest='MiC', default=25, help='The minimal zone coverage for which the second component of product giving the score will be maximal (integer), [default: %default]')
	parser.add_option( '', '--min_score', dest='min_score', default=70, help='The minimal score for a discordant zone to be identified as passed, [default: %default]')
	parser.add_option( '', '--out', dest='out', default='not_filled', help='Output file')
	parser.add_option( '', '--config', dest='config', default=None)
	(options, args) = parser.parse_args()

	pathname = os.path.dirname(sys.argv[0])

	loca_programs = ConfigParser.RawConfigParser()
	loca_programs.read(pathname+'/loca_programs.conf')

	if options.out == 'not_filled':
		mot = 'Please provide an argument for --out'
		raise ValueError(mot)
	if options.sam == 'not_filled':
		mot = 'Please provide an argument for --sam'
		raise ValueError(mot)

	logNameFile = os.path.splitext(options.sam)[0]+'.log'
	logOutput = open(logNameFile, 'w')
	if os.path.getsize(options.sam) == 0:
		outfile = open(options.out, 'w')
		outfile.close()
	else:
		tmp_name = tempfile.NamedTemporaryFile().name
		# tmp_name = 'toto'
		tmp_cov = tmp_name+'.cov'
		tmp_zone = tmp_name+'.zone'
		tmp_mate_zone = tmp_name+'_mate.zone'
		tmp_merge = tmp_name+'.merge'

		logOutput.write("\nsam file :"+str(options.sam))

		if options.config:
			config = ConfigParser.RawConfigParser()
			config.read(options.config)

			# Index bam file
			logOutput.write("\nstarting date : "+str(datetime.datetime.now()))
			t0 = datetime.datetime.now()
			index_bam_file(loca_programs, options.sam)
			logOutput.write("\nindex of the bam file : "+str(datetime.datetime.now() - t0))
			logOutput.flush()
			# Sort the bam upstream / downstream
			t0 = datetime.datetime.now()
			amont_aval(loca_programs, config.get('General','chr'), options.sam, config.get('Trie_discord','type'), tmp_name+'Sam')
			sam2bam(loca_programs, tmp_name+'Sam_amont', tmp_name+'Bam_amont')
			sam2bam(loca_programs, tmp_name+'Sam_aval', tmp_name+'Bam_aval')
			index_bam_file(loca_programs, tmp_name+'Bam_amont')
			index_bam_file(loca_programs, tmp_name+'Bam_aval')
			logOutput.write("\nSort and index amont/aval  : "+str(datetime.datetime.now() - t0))
			logOutput.flush()
			# Calcul the coverture site by site
			t0 = datetime.datetime.now()
			calcul_cov(loca_programs, options.sam, config.get('Trie_discord','type'), tmp_cov)
			logOutput.write('\ncalcul coverture : '+str(datetime.datetime.now() - t0))
			logOutput.flush()
			maxcov = config.getfloat('Calc_coverage','median_coverage')*config.getfloat('General','mult_max_cov')
			mincov = config.getfloat('Calc_coverage','median_coverage')*config.getfloat('General','mult_min_cov')
			minzone =  config.getint('General','min_zone')
			mingap = config.getint('General','min_gap')
			logOutput.write('\nMinimal accepted coverage:'+str(mincov))
			logOutput.write('\nMaximal accepted coverage:'+str(maxcov))
			logOutput.flush()

			# Find zones from the coverture
			t0 = datetime.datetime.now()
			statCovs = select_sur_couv(tmp_cov, config.getint('General','min_zone'), maxcov, mincov, config.getint('General','min_gap'), tmp_zone)
			nb_zone = statCovs[0]
			longerZone = statCovs[1]
			logOutput.write("\ntotal number of discordant zones : "+str(nb_zone))
			logOutput.write("\nlength of the longer zone : "+str(longerZone))
			logOutput.write('\nselect zones on coverture : '+str(datetime.datetime.now() - t0))
			logOutput.flush()
			os.remove(tmp_cov)
			ecart = config.getfloat('Calc_coverage','standard_deviation_insert')*3.0
			logOutput = open(logNameFile, 'a')
			logOutput.write('\nMargin:'+str(ecart))
			logOutput.write('\nNumber of zone to test:'+str(nb_zone))
			logOutput.flush()

			# Try to identify mate zones
			t0 = datetime.datetime.now()
			look_4_mate(logNameFile, loca_programs, config.get('Trie_discord','type'), options.sam, config.get('General','chr'), tmp_mate_zone, tmp_zone, ecart, config.getint('General','min_zone'), maxcov, mincov, mingap, nb_zone, longerZone, config.getfloat('Calc_coverage','median_insert'))
			logOutput.write('\nlook for mate zones : '+str(datetime.datetime.now() - t0))
			logOutput.flush()
			os.remove(tmp_zone)

			# Try to merge the mate zones identified above
			t0 = datetime.datetime.now()
			merge_zone(loca_programs, config.get('General','chr'), tmp_mate_zone, config.getfloat('General','max_dist_merge'), tmp_name+'Bam_amont', tmp_name+'Bam_aval', 'bam', tmp_merge)
			logOutput.write('\nmerge the zones : '+str(datetime.datetime.now() - t0))
			logOutput.flush()
			os.remove(tmp_mate_zone)

			# Calcul the score for each mate zones
			t0 = datetime.datetime.now()
			calculate_score(tmp_merge, config.getfloat('General','YiS'), config.getfloat('Score_discord','MiS'), config.getfloat('General','YiC'), config.getfloat('Score_discord','MiC'), config.getfloat('General','min_score'), config.get('General','chr'), options.out)
			logOutput.write('\ncalculate score : '+str(datetime.datetime.now() - t0))
			logOutput.flush()

			os.remove(tmp_merge)
			# Add header to the score file
			os.system("sed -i '1i#CHR-zone1\tSTART\tEND\tSIZE\tCOV\tCHR-zone2\tSTART\tEND\tSIZE\tCOV\tMISC\tREAD\tSCORE\tSTATUS' %s" % options.out)
		else:
			if options.ref == 'not_filled':
				raise ValueError('Please provide an argument for --ref')

			# Index bam file
			logOutput.write("\nstarting date : "+str(datetime.datetime.now()))
			t0 = datetime.datetime.now()
			index_bam_file(loca_programs, options.sam)
			logOutput.write("\nindex of the bam file : "+str(datetime.datetime.now() - t0))
			logOutput.flush()
			# Sort the bam upstream / downstream
			t0 = datetime.datetime.now()
			amont_aval(loca_programs, options.chr, options.sam, options.type, tmp_name+'Sam')
			sam2bam(loca_programs, tmp_name+'Sam_amont', tmp_name+'Bam_amont')
			sam2bam(loca_programs, tmp_name+'Sam_aval', tmp_name+'Bam_aval')
			index_bam_file(loca_programs, tmp_name+'Bam_amont')
			index_bam_file(loca_programs, tmp_name+'Bam_aval')
			logOutput.write("\nSort and index amont/aval  : "+str(datetime.datetime.now() - t0))
			logOutput.flush()

			# Calcul the coverture site by site
			t0 = datetime.datetime.now()
			calcul_cov(loca_programs, options.sam, options.type, tmp_cov)
			logOutput.write('\ncalcul coverture : '+str(datetime.datetime.now() - t0))
			logOutput.flush()
			maxcov = int(options.median_coverage)*int(options.mult_max_cov)
			mincov = int(options.median_coverage)*int(mult_min_cov)
			minzone =  options.min_zone
			mingap = options.min_gap
			logOutput.write('\nMinimal accepted coverage:'+str(mincov))
			logOutput.write('\nMaximal accepted coverage:'+str(maxcov))
			logOutput.flush()

			# Find zones from the coverture
			t0 = datetime.datetime.now()
			statCovs = select_sur_couv(tmp_cov, int(options.min_zone), float(options.maxcov), float(options.mincov), int(options.min_gap), tmp_zone)
			nb_zone = statCovs[0]
			longerZone = statCovs[1]
			logOutput.write("\ntotal number of discordant zones : "+str(nb_zone))
			logOutput.write("\nlength of the longer zone : "+str(longerZone))
			logOutput.write('\nselect zones on coverture : '+str(datetime.datetime.now() - t0))
			logOutput.flush()
			os.remove(tmp_cov)
			ecart = options.ecart
			logOutput.write('\nMargin:'+str(ecart))
			logOutput.write('\nNumber of zone to test:'+str(nb_zone))
			logOutput.flush()

			# Try to identify mate zones
			t0 = datetime.datetime.now()
			look_4_mate(logNameFile, loca_programs,  options.type, options.sam, options.chr, tmp_mate_zone, tmp_zone, float(options.ecart), int(options.min_zone), float(options.maxcov), float(options.mincov), int(options.min_gap), nb_zone, longerZone, int(options.med_insert))
			logOutput.write('\nlook for mate zones : '+str(datetime.datetime.now() - t0))
			logOutput.flush()
			os.remove(tmp_zone)

			# Try to merge the mate zones identified above
			t0 = datetime.datetime.now()
			merge_zone(loca_programs, options.chr, tmp_mate_zone, int(options.max_dist_merge), tmp_name+'Bam_amont', tmp_name+'Bam_aval', 'bam', tmp_merge)
			logOutput.write('\nmerge the zones : '+str(datetime.datetime.now() - t0))
			logOutput.flush()
			os.remove(tmp_mate_zone)

			# Calcul the score for each mate zones
			t0 = datetime.datetime.now()
			calculate_score(tmp_merge, float(options.YiS), float(options.MiS), float(options.YiC), float(optionsMiC), float(options.min_score), options.chr, options.out)
			logOutput.write('\ncalculate score : '+str(datetime.datetime.now() - t0))
			logOutput.flush()

			os.remove(tmp_merge)
			# Add header to the score file
			os.system("sed -i '1i#CHR-zone1\tSTART\tEND\tSIZE\tCOV\tCHR-zone2\tSTART\tEND\tSIZE\tCOV\tMISC\tREAD\tSCORE\tSTATUS' %s" % options.out)

	logOutput.write("\ntotal time : "+str(datetime.datetime.now() - t_start))
	logOutput.close()
	# if os.path.exists(options.sam+".bai"):
		# os.remove(options.sam+".bai")
	os.remove(logNameFile)
if __name__ == "__main__": __main__()
