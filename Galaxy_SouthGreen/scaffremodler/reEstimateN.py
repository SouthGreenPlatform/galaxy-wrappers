
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

import optparse, os, shutil, subprocess, sys, tempfile, fileinput, ConfigParser, operator, time, threading, datetime, multiprocessing

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def stop_err( msg ):
    sys.stderr.write( "%s\n" % msg )
    sys.exit()

def run_job (cmd_line, ERROR, ID):
	print cmd_line
	try:
		tmp = (tempfile.NamedTemporaryFile().name)+'-'+ID+'.error'
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


def recal_ins(FILE):
	LIST = []
	fichier = open(FILE)
	i = 0
	for line in fichier:
		data = line.split()
		if data != []:
			if line[0] != '@':
				if i == 0:
					i = 1
					LIST.append(abs(int(data[8])))
				else:
					i = 0
	if LIST == []:
		return 'NA'
	else:
		return sum(LIST)/len(LIST)
	
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


def __main__():
	#Parse Command Line
	parser = optparse.OptionParser(usage="python %prog [options]\n\nProgram designed by Guillaume MARTIN : guillaume.martin@cirad.fr\n\n"
	"This program re-estimate N present in DNA sequence. Re-estimated N are replaced by S")
	# Wrapper options.
	parser.add_option( '', '--config', dest='config', default='not_filled', help='A config file generated in ApMap pipeline')
	parser.add_option( '', '--exclude', dest='exclude', default='No_exclude', help='Chromosome names separated with "=" to exclude for the insert size estimation')
	parser.add_option( '', '--min_read', dest='min_read', default=30, help='The minimal read number requested to make the estimation, [default: %default]')
	parser.add_option( '', '--out', dest='out', default='N_restimated.fasta', help='The output name, [default: %default]')
	parser.add_option( '', '--thread', dest='thread', default=1, help='Number of thread to use, [default: %default]')
	(options, args) = parser.parse_args()
	
	
	
	ScriptPath = os.path.dirname(sys.argv[0])
	
	loca_programs = ConfigParser.RawConfigParser()
	loca_programs.read(ScriptPath+'/loca_programs.conf')
	
	if options.config == 'not_filled':
		sys.exit('--config argument is missing')
	if options.out == 'not_filled':
		sys.exit('--out argument is missing')
	
	config = ConfigParser.RawConfigParser()
	config.read(options.config)
	
	#identification of chromosomes to exclude
	CHR  = options.exclude.split('=')
	
	#Calculation of chromosomes size
	#1)Loading sequences
	record_dict = SeqIO.index(config.get('General','ref'), "fasta")
	#2)Recording informations in a temporary file
	f = tempfile.NamedTemporaryFile()
	liste_id = []
	for n in record_dict:
		if not(n in CHR):
			f.write(n+'\t1'+'\t'+str(len(str(record_dict[n].seq)))+'\n')
			liste_id.append(n)
			f.flush()
	#3)Generating a filtered sam from which insertsize will be estimated (on well mapped reads)
	tempo = tempfile.NamedTemporaryFile().name
	parse1 = '%s view -uf 2 -L %s %s | %s view -h -o %s - ' % (loca_programs.get('Programs','samtools'), f.name, config.get('Remove_dup','out'), loca_programs.get('Programs','samtools'), tempo+'_reEstimateN_mapped.sam')
	run_job(parse1, 'bug in parsing 1', "sam")
	f.close()
	INSERT = recal_ins(tempo+'_reEstimateN_mapped.sam')
	os.system('echo "Estimated insert size : '+str(INSERT)+'"')
	os.remove(tempo+'_reEstimateN_mapped.sam')
	# INSERT = 5400
	
	outseq = open(options.out,'w')
	#4)Estimation of the number of N for each N region for each chromosome
	# proc = int(os.popen('grep -c cores /proc/cpuinfo').read().split()[0])
	proc = int(options.thread)
	liste = []
	liste_process = []
	for n in liste_id:
		if not(n in CHR):
			outfile_seq = open(options.out+'_'+n+'_for_identSV.fasta','w')
			SeqIO.write(SeqRecord(record_dict[n].seq, id = n, description=''),outfile_seq, "fasta")
			outfile_seq.close()
			command_line = "%s %s/estimate.py --config %s --min_read %s --fasta %s --seq %s --out %s --insert %s" % (loca_programs.get('Programs','python'), ScriptPath, options.config, str(options.min_read), options.out+'_'+n+'_for_identSV.fasta', n, options.out+'_'+n+'_interm.fasta', str(INSERT))
			# t = multiprocessing.Process(target=restim, args=(tempo, config.get('Remove_dup','out'), config.get('General','orient'), int(options.min_read), str(record_dict[n].seq), n, options.out, INSERT,))
			t = multiprocessing.Process(target=run_job, args=(command_line, 'bug in estimate.py', n,))
			# Sticks the thread in a list so that it remains accessible 
			liste_process.append(t)
			liste.append(n)
			if len(liste) == proc:
				# Starts threads
				for process in liste_process:
					process.start()
				# This blocks the calling thread until the thread whose join() method is called is terminated.
				for process in liste_process:
					process.join()
				#the processes are done
				for k in liste:
					os.remove(options.out+'_'+k+'_for_identSV.fasta')
					record_fasta = SeqIO.index(options.out+'_'+k+'_interm.fasta', "fasta")
					SeqIO.write(SeqRecord(record_fasta[k].seq, id = k, description=''),outseq, "fasta")
					os.remove(options.out+'_'+k+'_interm.fasta')
					del record_fasta
				#remove tested reference 
				liste = []
				liste_process = []
	if liste:
		# Starts threads
		for process in liste_process:
			process.start()
		# This blocks the calling thread until the thread whose join() method is called is terminated.
		for process in liste_process:
			process.join()
		#the processes are done
		for k in liste:
			os.remove(options.out+'_'+k+'_for_identSV.fasta')
			record_fasta = SeqIO.index(options.out+'_'+k+'_interm.fasta', "fasta")
			SeqIO.write(SeqRecord(record_fasta[k].seq, id = k, description=''),outseq, "fasta")
			os.remove(options.out+'_'+k+'_interm.fasta')
			del record_fasta
		#remove tested reference 
		liste = []
		liste_process = []
	#5)Outputing sequence where N where not estimated
	for n in CHR:
		if n in record_dict:
			SeqIO.write(SeqRecord(record_dict[n].seq, id = k, description=''),outseq, "fasta")
	outseq.close()

if __name__ == "__main__": __main__()