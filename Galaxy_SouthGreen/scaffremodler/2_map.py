
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

import optparse, os, shutil, subprocess, sys, tempfile, fileinput, ConfigParser, operator, time, random, datetime, glob, sys
import multiprocessing

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


def worker(job):
	try:
		print "---"
		print job[0]
		print job[1]
		print "---"
		sys.stdout.flush()
		run_job(job[0], job[1])
	except Exception, e:
		print "error : "+e.__doc__+" ('"+e.message+")' in '"+job[0]+"'"



def Mapping(LOCA_PROGRAMS, TOOL, REF, Q1, Q2, ORIENT, MIN, MAX, QUAL, INDEX, RMINDEX, THREAD, OUT, PATHNAME):
		interm1 = OUT+'_mate1.sam'
		interm2 = OUT+'_mate2.sam'
		interm_sort1 = OUT+'_mate1_sorted.sam'
		interm_sort2 = OUT+'_mate2_sorted.sam'
		t0 = datetime.datetime.now()
		print t0
		if TOOL == 'bowtie':
			if INDEX == 'y':
				build_index = '%s -q %s %s' % (LOCA_PROGRAMS.get('Programs','bowtie-build'), REF, REF)
				run_job(build_index,'Indexing error:\n')
			if QUAL == '33':
				mapping1 = '%s --quiet -a -m 1 %s --phred33-quals -q %s -S %s' % (LOCA_PROGRAMS.get('Programs','bowtie'), REF, Q1, interm1)
				mapping2 = '%s --quiet -a -m 1 %s --phred33-quals -q %s -S %s' % (LOCA_PROGRAMS.get('Programs','bowtie'), REF, Q2, interm2)
				# mapping1 = '%s --quiet %s --phred33-quals -q %s -S %s' % (LOCA_PROGRAMS.get('Programs','bowtie'), REF, Q1, interm1)
				# mapping2 = '%s --quiet %s --phred33-quals -q %s -S %s' % (LOCA_PROGRAMS.get('Programs','bowtie'), REF, Q2, interm2)
			elif QUAL == '64':
				mapping1 = '%s --quiet -a -m 1 %s --solexa1.3-quals -q %s -S %s' % (LOCA_PROGRAMS.get('Programs','bowtie'), REF, Q1, interm1)
				mapping2 = '%s --quiet -a -m 1 %s --solexa1.3-quals -q %s -S %s' % (LOCA_PROGRAMS.get('Programs','bowtie'), REF, Q2, interm2)
				# mapping1 = '%s --quiet %s --solexa1.3-quals -q %s -S %s' % (LOCA_PROGRAMS.get('Programs','bowtie'), REF, Q1, interm1)
				# mapping2 = '%s --quiet %s --solexa1.3-quals -q %s -S %s' % (LOCA_PROGRAMS.get('Programs','bowtie'), REF, Q2, interm2)
			else:
				sys.exit('Unknown quality encoding : support only +33 or +64 encoding')
			sorting1 = '%s -jar %s SortSam INPUT=%s OUTPUT=%s SORT_ORDER=queryname QUIET=true MAX_RECORDS_IN_RAM=5000000 VERBOSITY=WARNING VALIDATION_STRINGENCY=SILENT' % (LOCA_PROGRAMS.get('Programs','java'), LOCA_PROGRAMS.get('Programs','picard-tool'), interm1, interm_sort1)
			sorting2 = '%s -jar %s SortSam INPUT=%s OUTPUT=%s SORT_ORDER=queryname QUIET=true MAX_RECORDS_IN_RAM=5000000 VERBOSITY=WARNING VALIDATION_STRINGENCY=SILENT' % (LOCA_PROGRAMS.get('Programs','java'), LOCA_PROGRAMS.get('Programs','picard-tool'), interm2, interm_sort2)
			merging = '%s %s/merge_sam.py --file1 %s --file2 %s --out %s --min %s --max %s --orient %s' % (LOCA_PROGRAMS.get('Programs','python'), PATHNAME, interm_sort1, interm_sort2, OUT, MIN, MAX, ORIENT)
			run_job(mapping1, 'Mapping error:\n')
			run_job(sorting1, 'Sorting error:\n')
			os.remove(interm1)
			run_job(mapping2, 'Mapping error:\n')
			run_job(sorting2, 'Sorting error:\n')
			os.remove(interm2)
			run_job(merging, 'Merging error:\n')
			os.remove(interm_sort1)
			os.remove(interm_sort2)
		elif TOOL == 'bowtie2_single':
			if INDEX == 'y':
				build_index = '%s -q -f %s %s' % (LOCA_PROGRAMS.get('Programs','bowtie2-build'), REF, REF)
				run_job(build_index,'Indexing error')
			mapping1 = '%s -D 20 -R 3 -N 0 -L 20 -i S,1,0.50 -x %s -q %s --phred%s -p %s -S %s' % (LOCA_PROGRAMS.get('Programs','bowtie2'), REF, Q1, QUAL, THREAD, interm1)
			mapping2 = '%s -D 20 -R 3 -N 0 -L 20 -i S,1,0.50 -x %s -q %s --phred%s -p %s -S %s' % (LOCA_PROGRAMS.get('Programs','bowtie2'), REF, Q2, QUAL, THREAD, interm2)
			sorting1 = '%s -jar %s SortSam INPUT=%s OUTPUT=%s SORT_ORDER=queryname QUIET=true MAX_RECORDS_IN_RAM=5000000 VERBOSITY=WARNING VALIDATION_STRINGENCY=SILENT' % (LOCA_PROGRAMS.get('Programs','java'), LOCA_PROGRAMS.get('Programs','picard-tool'), interm1, interm_sort1)
			sorting2 = '%s -jar %s SortSam INPUT=%s OUTPUT=%s SORT_ORDER=queryname QUIET=true MAX_RECORDS_IN_RAM=5000000 VERBOSITY=WARNING VALIDATION_STRINGENCY=SILENT' % (LOCA_PROGRAMS.get('Programs','java'), LOCA_PROGRAMS.get('Programs','picard-tool'), interm2, interm_sort2)
			merging = '%s %s/merge_sam.py --file1 %s --file2 %s --out %s --min %s --max %s --orient %s' % (LOCA_PROGRAMS.get('Programs','python'), PATHNAME, interm_sort1, interm_sort2, OUT, MIN, MAX, ORIENT)

			run_job(mapping1, 'Mapping error:\n')
			run_job(mapping2, 'Mapping error:\n')

			run_job(sorting1, 'Sorting error:\n')
			run_job(sorting2, 'Sorting error:\n')

			run_job(merging, 'Merging error:\n')
			os.remove(interm1)
			os.remove(interm2)
			os.remove(interm_sort1)
			os.remove(interm_sort2)
		elif TOOL == 'bowtie2':
			if INDEX == 'y':
				build_index = '%s -q -f %s %s' % (LOCA_PROGRAMS.get('Programs','bowtie2-build'), REF, REF)
				run_job(build_index,'Indexing error')
			mapping = '%s -D 20 -R 3 -N 0 -L 20 -i S,1,0.50 -I %s -X %s --%s -x %s -q1 %s -q2 %s --phred%s -p %s -S %s' % (LOCA_PROGRAMS.get('Programs','bowtie2'), MIN, MAX, ORIENT, REF, Q1, Q2, QUAL, THREAD, OUT)
			run_job(mapping,'Mapping error:\n')
		elif TOOL == 'bwa_mem':
			if INDEX == 'y':
				build_index = '%s index -a bwtsw %s 2>/dev/null' % (LOCA_PROGRAMS.get('Programs','bwa'), REF)
				run_job(build_index,'Indexing error:')
			mapping1 = '%s mem -t %s -M %s %s > %s' % (LOCA_PROGRAMS.get('Programs','bwa'), THREAD, REF, Q1, interm1)
			mapping2 = '%s mem -t %s -M %s %s > %s' % (LOCA_PROGRAMS.get('Programs','bwa'), THREAD, REF, Q2, interm2)
			sorting1 = '%s -jar %s SortSam INPUT=%s OUTPUT=%s SORT_ORDER=queryname QUIET=true MAX_RECORDS_IN_RAM=5000000 VERBOSITY=WARNING VALIDATION_STRINGENCY=SILENT' % (LOCA_PROGRAMS.get('Programs','java'), LOCA_PROGRAMS.get('Programs','picard-tool'), interm1, interm_sort1)
			sorting2 = '%s -jar %s SortSam INPUT=%s OUTPUT=%s SORT_ORDER=queryname QUIET=true MAX_RECORDS_IN_RAM=5000000 VERBOSITY=WARNING VALIDATION_STRINGENCY=SILENT' % (LOCA_PROGRAMS.get('Programs','java'), LOCA_PROGRAMS.get('Programs','picard-tool'), interm2, interm_sort2)
			merging = '%s %s/merge_sam.py --file1 %s --file2 %s --out %s --min %s --max %s --orient %s' % (LOCA_PROGRAMS.get('Programs','python'), PATHNAME, interm_sort1, interm_sort2, OUT, MIN, MAX, ORIENT)
			run_job(mapping1, 'Mapping error:\n')
			run_job(sorting1, 'Sorting error:\n')
			os.remove(interm1)
			run_job(mapping2, 'Mapping error:\n')
			run_job(sorting2, 'Sorting error:\n')
			os.remove(interm2)
			run_job(merging, 'Merging error:\n')
			os.remove(interm_sort1)
			os.remove(interm_sort2)
		elif TOOL == 'bwa':
			if INDEX == 'y':
				build_index = '%s index -a bwtsw %s 2>/dev/null' % (LOCA_PROGRAMS.get('Programs','bwa'), REF)
				run_job(build_index,'Indexing error:')
			bwasai1 = OUT+'_mate1.sai'
			bwasai2 = OUT+'_mate2.sai'
			mapping1 = '%s aln -t %s %s %s > %s' % (LOCA_PROGRAMS.get('Programs','bwa'), THREAD, REF, Q1, bwasai1)
			mapping2 = '%s aln -t %s %s %s > %s' % (LOCA_PROGRAMS.get('Programs','bwa'), THREAD, REF, Q2, bwasai2)
			mapping2_1 = '%s sampe %s %s %s %s %s > %s' % (LOCA_PROGRAMS.get('Programs','bwa'), REF, bwasai1, bwasai2, Q1, Q2, OUT)
			run_job(mapping1, 'Mapping error:\n')
			run_job(mapping2, 'Mapping error:\n')
			run_job(mapping2_1, 'Mapping2 error:\n')
			os.remove(bwasai1)
			os.remove(bwasai2)
		if RMINDEX == 'y':
			for filename in glob.glob(REF+'.*'):
				# print filename
				os.remove(filename)
		print datetime.datetime.now() - t0

def __main__():
	#Parse Command Line
	parser = optparse.OptionParser(usage="python %prog [options]\n\nProgram designed by Guillaume MARTIN : guillaume.martin@cirad.fr"
	"\n\n This script map paired reads on a reference and output a sam file containing the paired reads mapped.")
	# Wrapper options.
	parser.add_option( '', '--tool', dest='tool', default='bowtie2_single', help='The tool used : bowtie, bowtie2, bowtie2_single, bwa, bwa_mem, [default: %default]')
	parser.add_option( '', '--ref', dest='ref', default='not_filled', help='The multifasta reference file')
	parser.add_option( '', '--q1', dest='q1', default='not_filled', help='The mate1 fastq file')
	parser.add_option( '', '--q2', dest='q2', default='not_filled', help='The mate2 fastq file')
	parser.add_option( '', '--orient', dest='orient', default='rf', help='The expected orientation: rf or fr, [default: %default]')
	parser.add_option( '', '--mini', dest='mini', default='2500', help='The minimum insert size (integer), [default: %default]')
	parser.add_option( '', '--maxi', dest='maxi', default='7500', help='The maximum insert size (integer), [default: %default]')
	parser.add_option( '', '--qual', dest='qual', default='33', help='Fastq quality encoding: 33 or 64, [default: %default]')
	parser.add_option( '', '--index', dest='index', default='y', help='Build reference index : y or n,  [default: %default]')
	parser.add_option( '', '--rmindex', dest='rmindex', default='y', help='Remove reference index at the end of calculation: y or n, [default: %default]')
	parser.add_option( '', '--thread', dest='thread', default='1', help='The thread number used for mapping (integer), [default: %default]')
	parser.add_option( '', '--out', dest='out', default='mate.sam', help='The ouput of mapped reads, [default: %default]')
	parser.add_option( '', '--config', dest='config', default=None)
	(options, args) = parser.parse_args()



	pathname = os.path.dirname(sys.argv[0])

	loca_programs = ConfigParser.RawConfigParser()
	loca_programs.read(pathname+'/loca_programs.conf')

	if options.config:
		config = ConfigParser.RawConfigParser()
		config.read(options.config)
		Mapping(loca_programs, config.get('General','tool'), config.get('General','ref'), config.get('General','q1'), config.get('General','q2'), config.get('General','orient'), config.get('General','mini'), config.get('General','maxi'), config.get('General','qual'), config.get('General','index'), config.get('General','rmindex'), config.get('General','thread'), options.out, pathname)
		config.set('Mapping', 'out', options.out)
		if config.get('General','tool') in ['bowtie', 'bowtie2', 'bowtie2_single']:
			config.set('Single_filter', 'asxs', 1)
			config.set('Single_filter', 'qual', 'not_filled')
		else:
			config.set('Single_filter', 'asxs', 'not_filled')
			config.set('Single_filter', 'qual', 0)
		with open(options.config, 'wb') as configfile:
			config.write(configfile)
	else:#For commande line
		if options.q1 == 'not_filled':
			sys.exit('--q1 argument is missing')
		if options.q2 == 'not_filled':
			sys.exit('--q2 argument is missing')
		if options.ref == 'not_filled':
			sys.exit('--ref argument is missing')
		Mapping(loca_programs, options.tool, options.ref, options.q1, options.q2, options.orient, options.mini, options.maxi, options.qual, options.index, options.rmindex, options.thread, options.out, pathname)





if __name__ == "__main__": __main__()
