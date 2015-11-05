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

import optparse, os, shutil, subprocess, sys, tempfile, fileinput, ConfigParser, operator, time, multiprocessing

def stop_err( msg ):
    sys.stderr.write( "%s\n" % msg )
    sys.exit()

def run_job (cmd_line, ERROR):
	# print cmd_line
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



def __main__():
	#Parse Command Line
	parser = optparse.OptionParser()
	# Wrapper options. 
	parser.add_option( '', '--blast', dest='blast', default='not_filled', help='The blast file default output format')
	parser.add_option( '', '--ident', dest='ident', default='90', help='The minimal identity percentage, [default: %default]')
	parser.add_option( '', '--max_hit', dest='max_hit', default='1', help='The maximal hit number, [default: %default]')
	parser.add_option( '', '--seq', dest='seq', default='100', help='The number of sequence to search simultaneously, [default: %default]')
	parser.add_option( '', '--thread', dest='thread', default='1', help='The thread number used for mapping (integer), [default: %default]')
	parser.add_option( '', '--out', dest='out', default='not_filled', help='An id for intermediate output')
	(options, args) = parser.parse_args()

	ScriptPath = os.path.dirname(sys.argv[0])
	
	loca_programs = ConfigParser.RawConfigParser()
	loca_programs.read(ScriptPath+'/loca_programs.conf')
	
	proc = int(options.thread)
	
	if options.blast == 'not_filled':
		mot = 'Please provide an argument for --blast'
		sys.exit(mot)
	if options.out == 'not_filled':
		mot = 'Please provide an argument for --out'
		sys.exit(mot)
	
	file = open (options.blast)
	j = 0
	liste_job = []
	liste_temp = []
	for line in file:
		if line.split() != []:
			if line.split()[0] == 'BLASTN':
				if j == 0:
					outfile = open(options.out+'_File'+str(j)+'.temp','w')
					outfile.write(line)
					liste_job.append('%s %s/blat_results_analyzer_v3.pl --blat %s --identity %s --nombre %s > %s' % (loca_programs.get('Programs','perl'), ScriptPath, options.out+'_File'+str(j)+'.temp', options.ident, options.max_hit, options.out+'.bra'+str(j)))
					liste_temp.append(str(j))
				elif j%int(options.seq) == 0:
					outfile.close()
					liste_job.append('%s %s/blat_results_analyzer_v3.pl --blat %s --identity %s --nombre %s > %s' % (loca_programs.get('Programs','perl'), ScriptPath, options.out+'_File'+str(j)+'.temp', options.ident, options.max_hit, options.out+'.bra'+str(j)))
					liste_temp.append(str(j))
					outfile = open(options.out+'_File'+str(j)+'.temp','w')
					outfile.write(line)
				else:
					outfile.write(line)
				j += 1
			else:
				outfile.write(line)
		else:
			outfile.write(line)
	outfile.close()
	
	liste_process = []
	for n in liste_job:
		t = multiprocessing.Process(target=run_job, args=(n, 'Bug lauching blat_results_analyzer_v3.pl',))
		liste_process.append(t)
		if len(liste_process) == proc:
			# Starts threads
			for process in liste_process:
				process.start()
			# This blocks the calling thread until the thread whose join() method is called is terminated.
			for process in liste_process:
				process.join()
			#the processes are done
			liste_process = []
	if liste_process:
		# Starts threads
		for process in liste_process:
			process.start()
		# This blocks the calling thread until the thread whose join() method is called is terminated.
		for process in liste_process:
			process.join()
		#the processes are done
		liste_process = []

	os.system('cat '+options.out+'.bra* > '+options.out+'_final')
	for n in liste_temp:
		os.remove(options.out+'_File'+n+'.temp')
		os.remove(options.out+'.bra'+n)
	
	doublon = set()
	dico = {}
	file = open(options.out+'_final')
	for line in file:
		data = line.split()
		if data:
			if data[0] in dico:
				doublon.add(data[0])
			else:
				dico[data[0]] = line
	file.close()
	
	for n in dico:
		if not(n in doublon):
			print('\t'.join(dico[n].split()))
	os.remove(options.out+'_final')
	
if __name__ == "__main__": __main__()

