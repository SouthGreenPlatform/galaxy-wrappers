
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

def Filter(LOCA_PROGRAMS, SAM, TYPE, SORT, OUT):
	bamfile = OUT+'_filtered.bam'
	sortedbam = OUT+'_sorted.bam'
	rmdupbam = OUT+'_rmdup.bam'
	rmdupmetrics = OUT+'_rmdupmetrics.bam'
	if TYPE == 'sam':
		filter = '%s view -bS %s | %s view -uF 4 - | %s view -uF 8 - > %s' % (LOCA_PROGRAMS.get('Programs','samtools'), SAM, LOCA_PROGRAMS.get('Programs','samtools'), LOCA_PROGRAMS.get('Programs','samtools'), bamfile)
	elif TYPE == 'bam':
		filter = '%s view -uF 4 %s | %s view -uF 8 - > %s' % (LOCA_PROGRAMS.get('Programs','samtools'), SAM, LOCA_PROGRAMS.get('Programs','samtools'), bamfile)
	else:
		mot = SAM+' argument passed in --sam is not recognized'
		sys.exit(mot)
	sorting1 = '%s -jar %s SortSam INPUT=%s OUTPUT=%s SORT_ORDER=coordinate QUIET=true MAX_RECORDS_IN_RAM=5000000 VERBOSITY=WARNING VALIDATION_STRINGENCY=SILENT' % (LOCA_PROGRAMS.get('Programs','java'), LOCA_PROGRAMS.get('Programs','picard-tool'), bamfile, sortedbam)
	rmdup = '%s -jar %s MarkDuplicates INPUT=%s OUTPUT=%s METRICS_FILE=%s REMOVE_DUPLICATES=true QUIET=true MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 VERBOSITY=WARNING VALIDATION_STRINGENCY=SILENT' % (LOCA_PROGRAMS.get('Programs','java'), LOCA_PROGRAMS.get('Programs','picard-tool'), sortedbam, rmdupbam, rmdupmetrics)
	sorting2 = '%s -jar %s SortSam INPUT=%s OUTPUT=%s SORT_ORDER=%s QUIET=true MAX_RECORDS_IN_RAM=5000000 VERBOSITY=WARNING VALIDATION_STRINGENCY=SILENT' % (LOCA_PROGRAMS.get('Programs','java'), LOCA_PROGRAMS.get('Programs','picard-tool'), rmdupbam, OUT, SORT)
	
	run_job(filter, 'Error in first filter:')
	run_job(sorting1, 'Error in sorting1:')
	os.remove(bamfile)
	run_job(rmdup, 'Error in removing duplicates:')
	os.remove(sortedbam)
	run_job(sorting2, 'Error in sorting2:')
	os.remove(rmdupbam)
	os.system("sed -n 8p %s | cut -f 3 | sed 's/^/Read pairs examined : /'" % rmdupmetrics)
	os.system("sed -n 8p %s | cut -f 6 | sed 's/^/Read pairs duplicates : /'" % rmdupmetrics)
	os.system("sed -n 8p %s | cut -f 7 | sed 's/^/Read pairs optical duplicates : /'" % rmdupmetrics)
	os.system("sed -n 8p %s | cut -f 8 | sed 's/^/Duplication proportion : /'" % rmdupmetrics)
	os.remove(rmdupmetrics)

def __main__():
	#Parse Command Line
	parser = optparse.OptionParser(usage="python %prog [options]\n\nProgram designed by Guillaume MARTIN : guillaume.martin@cirad.fr"
	"\n\n This script remove unmapped paired reads and read duplicates")
	# Wrapper options. 
	parser.add_option( '', '--sam', dest='sam', default='not_filled', help='Paired sam/bam file')
	parser.add_option( '', '--type', dest='type', default='sam', help='Input type : sam or bam, [default: %default]')
	parser.add_option( '', '--sort', dest='sort', default='coordinate', help='Sort order queryname or coordinate, [default: %default]')
	parser.add_option( '', '--rminput', dest='rminput', default='n', help='Remove input file: y or n, [default: %default]')
	parser.add_option( '', '--out', dest='out', default='rmdup_mapped.bam', help='Output file')
	parser.add_option( '', '--config', dest='config', default=None)
	(options, args) = parser.parse_args()
	
	
	
	pathname = os.path.dirname(sys.argv[0])
	
	loca_programs = ConfigParser.RawConfigParser()
	loca_programs.read(pathname+'/loca_programs.conf')
	
	if options.config:
		config = ConfigParser.RawConfigParser()
		config.read(options.config)
		if options.sam == 'not_filled':
			Filter(loca_programs, config.get('Single_filter','out'), config.get('Single_filter','type'), config.get('Remove_dup','sort'), options.out)
		else:
			Filter(loca_programs, options.sam, config.get('Single_filter','type'), config.get('Remove_dup','sort'), options.out)
		if config.get('Remove_dup','rminput') == 'y':
			os.remove(config.get('Single_filter','out'))
		config.set('Remove_dup', 'out', options.out)
		config.set('Remove_dup', 'type', 'bam')
		with open(options.config, 'wb') as configfile:
			config.write(configfile)
	else:
		if options.sam == 'not_filled':
			sys.exit('--sam argument is missing, please provide a bam or a sam')
		Filter(loca_programs, options.sam, options.type, options.sort, options.out)
		if options.rminput == 'y':
			os.remove(options.sam)
	
	
if __name__ == "__main__": __main__()



