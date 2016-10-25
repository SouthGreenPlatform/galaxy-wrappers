
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

import optparse, os, shutil, subprocess, sys, tempfile, fileinput, ConfigParser, operator, time, datetime

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

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
	
def estimateN(chromosome, debut, fin, OR, sam, debut_rec, fin_rec, MIN_READ):
	LIST = []
	fichier = open(sam)
	for line in fichier:
		data = line.split()
		if data != []:
			if line[0] != '@':
				if data[6] == '=' and data[2] == chromosome:
					if int(data[3]) <= int(data[7]):#for selection of read spanning the zone
						if int(data[3]) <= debut and fin <= int(data[7]) and debut_rec <= int(data[3]) and int(data[7]) <= fin_rec:
							if data[1] == '83' or data[1] == '163' or data[1] == '99' or data[1] == '147':
								#concordant reads
								LIST.append(abs(int(data[8])))
								# print data[0:9]
							elif data[1] == '81':#mate1 R et mate2 F et il s'agit de mate1
								if int(data[3]) < int(data[7]) and OR == 'rf':
									LIST.append(abs(int(data[8])))
									# print data[0:9]
								elif int(data[3]) > int(data[7]) and OR == 'fr':
									LIST.append(abs(int(data[8])))
									print data[0:9]
									sys.exit('There is a probleme in estimateN 1')
									# print data[0:9]
							elif data[1] == '161':#mate1 R et mate2 F et il s'agit de mate2
								if int(data[3]) > int(data[7]) and OR == 'rf':
									LIST.append(abs(int(data[8])))
									print data[0:9]
									sys.exit('There is a probleme in estimateN 2')
									# print data[0:9]
								elif int(data[3]) < int(data[7]) and OR == 'fr':
									LIST.append(abs(int(data[8])))
									# print data[0:9]
							elif data[1] == '97':#mate1 F et mate2 R et il s'agit de mate1
								if int(data[3]) > int(data[7]) and OR == 'rf':
									LIST.append(abs(int(data[8])))
									print data[0:9]
									sys.exit('There is a probleme in estimateN 3')
									# print data[0:9]
								elif int(data[3]) < int(data[7]) and OR == 'fr':
									LIST.append(abs(int(data[8])))
									# print data[0:9]
							elif data[1] == '145':#mate1 F et mate2 R et il s'agit de mate2
								if int(data[3]) < int(data[7]) and OR == 'rf':
									LIST.append(abs(int(data[8])))
									# print data[0:9]
								elif int(data[3]) > int(data[7]) and OR == 'fr':
									LIST.append(abs(int(data[8])))
									print data[0:9]
									sys.exit('There is a probleme in estimateN 4')
									# print data[0:9]
					else:
						if int(data[7]) <= debut and fin <= int(data[3]) and debut_rec <= int(data[7]) and int(data[3]) <= fin_rec:
							if data[1] == '83' or data[1] == '163' or data[1] == '99' or data[1] == '147':
								LIST.append(abs(int(data[8])))
								# print data[0:9]
							elif data[1] == '81':#mate1 R et mate2 F et il s'agit de mate1
								if int(data[3]) < int(data[7]) and OR == 'rf':
									LIST.append(abs(int(data[8])))
									print data[0:9]
									sys.exit('There is a probleme in estimateN 5')
									# print data[0:9]
								elif int(data[3]) > int(data[7]) and OR == 'fr':
									LIST.append(abs(int(data[8])))
									# print data[0:9]
							elif data[1] == '161':#mate1 R et mate2 F et il s'agit de mate2
								if int(data[3]) > int(data[7]) and OR == 'rf':
									LIST.append(abs(int(data[8])))
									# print data[0:9]
								elif int(data[3]) < int(data[7]) and OR == 'fr':
									LIST.append(abs(int(data[8])))
									print data[0:9]
									sys.exit('There is a probleme in estimateN 6')
									# print data[0:9]
							elif data[1] == '97':#mate1 F et mate2 R et il s'agit de mate1
								if int(data[3]) > int(data[7]) and OR == 'rf':
									LIST.append(abs(int(data[8])))
									# print data[0:9]
								elif int(data[3]) < int(data[7]) and OR == 'fr':
									LIST.append(abs(int(data[8])))
									print data[0:9]
									sys.exit('There is a probleme in estimateN 7')
									# print data[0:9]
							elif data[1] == '145':#mate1 F et mate2 R et il s'agit de mate2
								if int(data[3]) < int(data[7]) and OR == 'rf':
									LIST.append(abs(int(data[8])))
									print data[0:9]
									sys.exit('There is a probleme in estimateN 8')
									# print data[0:9]
								elif int(data[3]) > int(data[7]) and OR == 'fr':
									LIST.append(abs(int(data[8])))
									# print data[0:9]
	print len(LIST)
	if LIST == []:
		return 'NA'
	elif len(LIST) < (2*int(MIN_READ)):
		return 'NA'
	else:
		return sum(LIST)/len(LIST)

def limit(DEBUT, FIN, SEQUENCE):
	debut_recal = 0
	fin_recal = len(SEQUENCE)-1
	j = 0
	for k in SEQUENCE:
		if k == 'n' or k == 'N':
			if j < DEBUT:
				debut_recal = j+1
			elif j > FIN:
				# print j, '>', FIN
				fin_recal = j-1
				break
		j = j + 1
	return [debut_recal, fin_recal]

def restim(LOCA_PROGRAMS, TEMP, BAM, ORIENT, MIN_READ, SEQ, NAME, OUT, INSERT):
	SEQUENCE = ''
	#Par chromosome on extrait les reads mappant sur ce chromosome
	tempo = tempfile.NamedTemporaryFile()
	tempo.write(NAME+'\t1'+'\t'+str(len(str(SEQ)))+'\n')
	tempo.flush()
	parse = '%s view -b -L %s %s > %s' % (LOCA_PROGRAMS.get('Programs','samtools'), tempo.name, BAM, TEMP+'chr'+NAME+'_reEstimateN_mapped.bam')
	run_job(parse, 'bug in parsing')
	tempo.close()
	#On enregistre la sequence du chromosome
	sequence = str(SEQ)
	i = len(sequence)
	avant = i
	if sequence[avant-1] == 'n' or sequence[avant-1] == 'N':
		sys.exit('There is a bug in the script')
	debut = ''
	fin = ''
	#on travail en partant de la fin de la sequence
	# t0 = datetime.datetime.now()
	while i > 0:
		caractere = sequence[i-1:i]
		#1)On cherche une region contenant des N
		if caractere == 'n' or caractere == 'N':
			if fin == '':
				#on arrive sur une region avec des N
				fin = i - 1
				debut = i - 1
			else:
				#on est sur des N
				debut = i - 1
		elif debut != '' and fin != '':
			# t0 = datetime.datetime.now()
			SEQUENCE = sequence[fin+1:avant] + SEQUENCE
			avant = debut
			# print 'new_seq', datetime.datetime.now() - t0
			# t0 = datetime.datetime.now()
			#2)A piece of verification
			if sequence[debut] != 'n' and sequence[debut] != 'N':
				sys.exit('There is a bug in the script')
			if sequence[fin] != 'n' and sequence[fin] != 'N':
				sys.exit('There is a bug in the script')
			if sequence[debut-1] == 'n' or sequence[debut-1] == 'N':
				sys.exit('There is a bug in the script')
			if sequence[fin+1] == 'n' or sequence[fin+1] == 'N':
				sys.exit('There is a bug in the script')
			#end of the piece of verification
			#3)recalculation of 'debut' and 'fin'
			if (debut - (2*INSERT)) >= 0 and (fin+(2*INSERT)) <= (len(sequence) - 1):#no problem
				# print (2*INSERT), (2*INSERT)+fin-debut
				RECAL = limit((2*INSERT), ((2*INSERT)+fin-debut), sequence[debut-(2*INSERT):fin+(2*INSERT)+1])
				debut_recal = debut - (2*INSERT) + RECAL[0]
				fin_recal = debut - (2*INSERT) + RECAL[1]
			elif (debut - (2*INSERT)) >= 0 and (fin+(2*INSERT)) > (len(sequence) - 1):#with 2*INSERT we are over sequence length
				RECAL = limit((2*INSERT), ((2*INSERT)+fin-debut), sequence[debut-(2*INSERT):fin+(2*INSERT)+1])
				debut_recal = debut - (2*INSERT) + RECAL[0]
				fin_recal = debut - (2*INSERT) + RECAL[1]
			elif (debut - (2*INSERT)) < 0 and (fin+(2*INSERT)) <= (len(sequence) - 1):#with 2*INSERT we are in negative coordinate
				RECAL = limit(debut, fin, sequence[0:fin+(2*INSERT)+1])
				debut_recal = RECAL[0]
				fin_recal = RECAL[1]
			else:#with 2*INSERT first coordinate is negative and the second is over sequence length
				RECAL = limit(debut, fin, sequence)
				debut_recal = RECAL[0]
				fin_recal = RECAL[1]
			#4)On recupere les reads qui mappent dans la region
			f = tempfile.NamedTemporaryFile()
			if (debut - (2*INSERT)) >= 0 and (fin+(2*INSERT)) <= (len(sequence) - 1):#no problem
				f.write(NAME+'\t'+str((debut+1)-(2*INSERT))+'\t'+str((fin+1)+(2*INSERT))+'\n')
			elif (debut - (2*INSERT)) >= 0 and (fin+(2*INSERT)) > (len(sequence) - 1):#with 2*INSERT we are over sequence length
				f.write(NAME+'\t'+str((debut+1)-(2*INSERT))+'\t'+str(len(sequence))+'\n')
			elif (debut - (2*INSERT)) < 0 and (fin+(2*INSERT)) <= (len(sequence) - 1):#with 2*INSERT we are in negative coordinate
				f.write(NAME+'\t1\t'+str((fin+1)+(2*INSERT))+'\n')
			else:#with 2*INSERT first coordinate is negative and the second is over sequence length
				f.write(NAME+'\t1\t'+str(len(sequence))+'\n')
			f.flush()
			sub_sam = '%s view -u -L %s %s | %s view -h -o %s - ' % (LOCA_PROGRAMS.get('Programs','samtools'), f.name, TEMP+'chr'+NAME+'_reEstimateN_mapped.bam', LOCA_PROGRAMS.get('Programs','samtools'), TEMP+'chr'+NAME+'_reEstimateN_mapped.sam')
			run_job(sub_sam, 'bug in read filtering around zone')
			f.close()
			INSERT_est = estimateN(NAME, debut, fin, ORIENT, TEMP+'chr'+NAME+'_reEstimateN_mapped.sam', debut_recal, fin_recal, MIN_READ)
			if INSERT_est == 'NA':
				#on peut pas estimer les N
				# print NAME, debut, fin, fin-debut+1 ,INSERT, INSERT_est, 'NA', "(",debut_recal, fin_recal,")"
				print (NAME+'\t'+str(debut)+'\t'+str(fin)+'\t'+str(fin-debut+1)+'\t'+str(INSERT)+'\t'+str(INSERT_est)+'\tNA\n')
				#to put N in the sequence
				iter = 0
				while iter < ((fin - debut) + 1):
					SEQUENCE = 'N' + SEQUENCE
					iter = iter + 1
			else:
				#On peut estimer les N
				# print NAME, debut, fin, fin-debut+1 ,INSERT, INSERT_est, (fin-debut+1)-(INSERT_est - INSERT), "(",debut_recal, fin_recal,")"
				print (NAME+'\t'+str(debut)+'\t'+str(fin)+'\t'+str(fin-debut+1)+'\t'+str(INSERT)+'\t'+str(INSERT_est)+'\t'+str((fin-debut+1)-(INSERT_est - INSERT))+'\n')
				#to put N in the sequence
				if (fin-debut+1)-(INSERT_est - INSERT) > 0:
					#La region contient des N
					iter = 0
					while iter < (fin-debut+1)-(INSERT_est - INSERT):
						SEQUENCE = 'E' + SEQUENCE
						iter = iter + 1
				else:
					#On estime que les deux regions semblent plus proches que possible on met par defaut 20 S
					iter = 0
					while iter < 20:
						SEQUENCE = 'E' + SEQUENCE
						iter = iter + 1
			os.remove(TEMP+'chr'+NAME+'_reEstimateN_mapped.sam')
			debut = ''
			fin = ''
			# SEQUENCE = caractere + SEQUENCE
			# print 'manage N', datetime.datetime.now() - t0
		# else:
			# SEQUENCE = caractere + SEQUENCE
		i = i - 1
	SEQUENCE = sequence[0:avant] + SEQUENCE
	os.remove(TEMP+'chr'+NAME+'_reEstimateN_mapped.bam')
	outfile = open(OUT,'w')
	SeqIO.write(SeqRecord(Seq(SEQUENCE, generic_dna), id = NAME, description=''),outfile, "fasta")
	outfile.close()


def __main__():
	#Parse Command Line
	parser = optparse.OptionParser(usage="python %prog [options]\n\nProgram designed by Guillaume MARTIN : guillaume.martin@cirad.fr\n\n"
	"This program re-estimate N present in DNA sequence. Restimated N are replaced by S")
	# Wrapper options.
	parser.add_option( '', '--config', dest='config', default='not_filled', help='A config file generated in ApMap pipeline')
	parser.add_option( '', '--min_read', dest='min_read', default='not_filled', help='The minimal read number requested to make the estimation, [default: %default]')
	parser.add_option( '', '--fasta', dest='fasta', default='not_filled', help='The fasta file, [default: %default]')
	parser.add_option( '', '--seq', dest='seq', default='not_filled', help='The sequence name, [default: %default]')
	parser.add_option( '', '--out', dest='out', default='not_filled', help='The output name, [default: %default]')
	parser.add_option( '', '--insert', dest='insert', default='not_filled', help='Re-estimated insert size, [default: %default]')
	(options, args) = parser.parse_args()
	
	
	
	pathname = os.path.dirname(sys.argv[0])
	
	loca_programs = ConfigParser.RawConfigParser()
	loca_programs.read(pathname+'/loca_programs.conf')
	
	if options.config == 'not_filled':
		mot = 'Please provide an argument for --config'
		sys.exit(mot)
	if options.min_read == 'not_filled':
		mot = 'Please provide an argument for --min_read'
		sys.exit(mot)
	if options.fasta == 'not_filled':
		mot = 'Please provide an argument for --fasta'
		sys.exit(mot)
	if options.seq == 'not_filled':
		mot = 'Please provide an argument for --seq'
		sys.exit(mot)
	if options.out == 'not_filled':
		mot = 'Please provide an argument for --out'
		sys.exit(mot)
	if options.insert == 'not_filled':
		mot = 'Please provide an argument for --insert'
		sys.exit(mot)
	
	t0 = datetime.datetime.now()
	config = ConfigParser.RawConfigParser()
	config.read(options.config)

	record_dict = SeqIO.index(options.fasta, "fasta")
	
	tempo = tempfile.NamedTemporaryFile().name
	
	restim(loca_programs, tempo, config.get('Remove_dup','out'), config.get('General','orient'), int(options.min_read), str(record_dict[options.seq].seq), options.seq, options.out, int(options.insert))
	print 'done', datetime.datetime.now() - t0
if __name__ == "__main__": __main__()