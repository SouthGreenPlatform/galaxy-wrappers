
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

def recal_ins(LOCA_PROGRAMS, FILE, TYPE):
	# os.system('echo "Calculating insert size on the firts approximation of the well mapped read"')
	tmp = tempfile.NamedTemporaryFile().name
	# print tmp
	if TYPE == 'bam':
		estim_insert = '%s view -uf 2 %s | %s view -h -o %s - ' % (LOCA_PROGRAMS.get('Programs','samtools'), FILE, LOCA_PROGRAMS.get('Programs','samtools'), tmp)
	elif TYPE == 'sam':
		estim_insert = '%s view -bS %s | %s view -uf 2 - | %s view -h -o %s - ' % (LOCA_PROGRAMS.get('Programs','samtools'), FILE, LOCA_PROGRAMS.get('Programs','samtools'), LOCA_PROGRAMS.get('Programs','samtools'), tmp)
	else:
		mot = TYPE+' argument passed in --sam is not recognized'
		sys.exit(mot)
	
	run_job(estim_insert, 'Error in bam/sam filtering:\n')
	
	LIST = []
	fichier = open(tmp)
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
	MED = mediane(LIST)
	MOY = moyenne(LIST)
	EC = ecart_type(LIST)
	LIST = []
	os.remove(tmp)
	return [MED, MOY, EC]

#Fonction that calculate la median, mean and interval containing (INT*100)% of values, COL: column to treate (0 based)
def stat(FILE, COL, INT, VERBOSE, STAT_FILE):
	FICH = open(FILE)
	DIC = {}
	DIC_final = {}
	for LINE in FICH:
		DATA = LINE.split()
		if DATA != []:
			if DATA[0] in DIC:
				DIC[DATA[0]].append(float(DATA[COL]))
			else:
				DIC[DATA[0]] = []
				DIC[DATA[0]].append(float(DATA[COL]))
	for n in DIC:
		DIC_final[n] = [moyenne(DIC[n]), mediane(DIC[n]), intervalle(DIC[n], INT), len(DIC[n])]
	somme_moy = 0
	somme_med = 0
	taille = 0
	outfile = open(STAT_FILE, 'w')
	for n in DIC_final:
		if VERBOSE == 'all':
			outfile.write("Chromosome : "+n+" :\n")
			outfile.write("\tMean coverage             : "+str(DIC_final[n][0])+"\n")
			outfile.write("\tMedian coverage           : "+str(DIC_final[n][1])+"\n")
			outfile.write("\tConfidence interval (90%) : ["+str(DIC_final[n][2][0])+","+str(DIC_final[n][2][1])+"]\n")
		somme_moy = somme_moy + DIC_final[n][0]*DIC_final[n][3]
		somme_med = somme_med + DIC_final[n][1]*DIC_final[n][3]
		taille = taille + DIC_final[n][3]
	if VERBOSE == 'all' or VERBOSE == 'yes':
		outfile.write("*****Average*****\n")
		outfile.write("\tMean coverage   : "+str(somme_moy/taille)+"\n")
		outfile.write("\tMedian coverage : "+str(somme_med/taille)+"\n")
	outfile.close()
	return [somme_moy/taille, somme_med/taille]

def moyenne(L):
	if len(L) == 0:
		moyenne = 0
	else:
		moyenne = sum(L)/float(len(L))
	return moyenne

def mediane(L):
	L.sort()
	N = len(L)
	n = N/2.0
	p = int(n)
	if n == 0:
		return 0
	if n == 1:
		return (L[0])
	elif n == p:
		return (L[p-1]+L[p])/2.0
	else:
		return L[p]
		
def intervalle(L, P):
	L.sort()
	N = len(L)
	return [L[int((N-(N*P))/2.0)],L[int(N-((N-(N*P))/2.0))]]

def variance(L) :
	n = len(L)
	mq = moyenne(L)**2
	s = sum([x**2 for x in L])
	variance = (n/(n-1))*(s / n - mq)
	return variance
	
def ecart_type(L) :
	VAR = variance(L)
	ecart_type = math.sqrt(VAR)
	return ecart_type

def calcul_stat(LOCA_PROGRAMS, SAM, TYPE, OUT, STAT_FILE):
	if TYPE == 'bam':
		cal_cov = '%s depth %s > %s' % (LOCA_PROGRAMS.get('Programs','samtools'), SAM, OUT)
	elif TYPE == 'sam':
		sam2bam = '%s view -bS %s > %s' % (LOCA_PROGRAMS.get('Programs','samtools'), SAM, OUT+'.bam')
		run_job(sam2bam, 'Error in sam2bam conversion:\n')
		cal_cov = '%s depth %s > %s' %  (LOCA_PROGRAMS.get('Programs','samtools'), OUT+'.bam', OUT)
	else:
		mot = TYPE+' argument passed in --type is not recognized'
		sys.exit(mot)
	run_job(cal_cov, 'Error in calculating coverage:\n')
	if TYPE == 'sam':
		os.remove(OUT+'.bam')
	
	liste = stat(OUT, 2, 0.9, 'all', STAT_FILE)
	
	INFO_INSERT = recal_ins(LOCA_PROGRAMS, SAM, TYPE)
	insert = float(INFO_INSERT[0])
	standev = str(INFO_INSERT[2])
	outfile = open(STAT_FILE, 'a')
	outfile.write("Insert size as been re-estimated to :"+str(insert)+"\n")
	outfile.write("Standard deviation of insert size as been re-estimated to :"+standev+"\n")
	return [INFO_INSERT[0], INFO_INSERT[1], INFO_INSERT[2], liste[0], liste[1]]

	
def __main__():
	#Parse Command Line
	parser = optparse.OptionParser(usage="python %prog [options]\n\nProgram designed by Guillaume MARTIN : guillaume.martin@cirad.fr\n\nThis program take in input a sam/bam file,"
	" calculate coverage for each covered sites of the reference sequences, estimate mean, median and 90% confidence interval coverage for the covered sites and estimate mean and "
	"standard deviation of library insert size.")
	
	parser.add_option( '', '--sam', dest='sam', default='not_filled', help='Paired sam/bam file')
	parser.add_option( '', '--type', dest='type', default='sam', help='Input type : sam or bam, [default: %default]')
	parser.add_option( '', '--out', dest='out', default='coverage.cov', help='Output file')
	parser.add_option( '', '--stat', dest='stat', default='stat.txt', help='Output statistic file')
	parser.add_option( '', '--outconf', dest='outconf', default='stat.conf', help='Output configuration file with statistics')
	parser.add_option( '', '--config', dest='config', default=None)
	(options, args) = parser.parse_args()
	
	
	
	pathname = os.path.dirname(sys.argv[0])
	
	loca_programs = ConfigParser.RawConfigParser()
	loca_programs.read(pathname+'/loca_programs.conf')
	
	if options.config:
		config = ConfigParser.RawConfigParser()
		config.read(options.config)
		if options.sam == 'not_filled':
			STAT = calcul_stat(loca_programs, config.get('Remove_dup','out'), config.get('Remove_dup','type'), options.out, options.stat)
		else:
			STAT = calcul_stat(loca_programs, options.sam, config.get('Remove_dup','type'), options.out, options.stat)
		mini = float(STAT[0] - (config.getfloat('General','sd_multiplicator')*STAT[2]))
		maxi = float(STAT[0] + (config.getfloat('General','sd_multiplicator')*STAT[2]))
		config.set('Calc_coverage', 'out', options.out)
		config.set('Calc_coverage', 'median_insert', STAT[0])
		config.set('Calc_coverage', 'mean_insert', STAT[1])
		config.set('Calc_coverage', 'standard_deviation_insert', STAT[2])
		config.set('Calc_coverage', 'mean_coverage', STAT[3])
		config.set('Calc_coverage', 'median_coverage', STAT[4])
		config.set('Calc_coverage', 'mini', mini)
		config.set('Calc_coverage', 'maxi', maxi)
		config.set('Score_discord', 'MiS', (config.getfloat('General','MiS')*STAT[0]))
		config.set('Score_discord', 'MiC', (config.getfloat('General','MiC')*STAT[4]))
		with open(options.config, 'wb') as configfile:
			config.write(configfile)
		copy = 'cp %s %s' % (options.config, options.outconf)
		os.system(copy)
	else:
		if options.sam == 'not_filled':
			mot = 'Please provide an argument for --sam'
			sys.exit(mot)
		calcul_stat(loca_programs, options.sam, options.type, options.out, options.stat)
	
	
if __name__ == "__main__": __main__()