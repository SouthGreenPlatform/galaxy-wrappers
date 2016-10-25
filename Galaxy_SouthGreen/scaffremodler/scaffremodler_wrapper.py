
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

import optparse, os, shutil, subprocess, sys, tempfile, fileinput, ConfigParser, operator, time, random, multiprocessing, datetime, math, re


def stop_err( msg ):
	sys.stderr.write( "%s\n" % msg )
	sys.exit()

def run_job (cmd_line, ERROR):
	print str(datetime.datetime.now())+" : "+cmd_line
	sys.stdout.flush()
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
		os.system('rm '+tmp)
		if returncode != 0:
			raise Exception, stderr
	except Exception, e:
		stop_err( ERROR + str( e ) )
	finally:
		return returncode


def main(job):

	try:
		#Modify the different paths for working in sub folders.
		regexSamFile = re.compile("(--sam\s.+)")
		regexOutFile = re.compile("(--out\s.+)")
		samFile = os.path.splitext(regexSamFile.search(job).group(1).split()[1])[0]
		outFile = regexOutFile.search(job).group(1).split()[1]
		job = job.replace('--sam ','--sam ../')
		job = job.replace('--config ','--config ../')

		# Create a dir based on the discordant type
		WORKING_DIR = tempfile.mkdtemp(prefix="temp_"+samFile+'_', dir=os.getcwd()).split('/')[-1]+'/'
		os.chdir(WORKING_DIR)
		returnCode = run_job(job, "error on the job : "+job+"\n")
	except Exception as e:
		print e
	finally:
		if os.path.isfile(outFile):
			os.rename(outFile, os.pardir+"/"+outFile)
		os.chdir(os.pardir)
		shutil.rmtree(WORKING_DIR)
		return returnCode


# def worker(listJobs, out_q):
	# """

	# """
	# outdict = {}
	# # print "###"
	# # print listJobs
	# # print "###"
	# for job in listJobs:
		# try:
			# # print "job : "+job
			# outdict[job] = main(job)
		# except Exception as e:
			# outdict[job] = 0
			# sys.stderr.write(format(e))
			# pass
	# out_q.put(outdict)

def __main__():
	#Parse Command Line
	parser = optparse.OptionParser(usage="python %prog [options]\n\nProgram designed by Guillaume MARTIN : guillaume.martin@cirad.fr\n\n"
	"Is a wrapper for ApMap tools")
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
	parser.add_option( '', '--filter_multi', dest='filter_multi', default='y', help='Filter reads with multiple locations : y or n,  [default: %default]')
	parser.add_option( '', '--mini_dis', dest='mini_dis', default='10000', help='The minimal insert size to keep the discordant read for structural variation search (integer), [default: %default]')
	parser.add_option( '', '--mult_max_cov', dest='mult_max_cov', default='10', help='multiplicator of median coverage for maximal median coverage to keep a zone (float), [default: %default]')
	parser.add_option( '', '--mult_min_cov', dest='mult_min_cov', default='0.25', help='multiplicator of median coverage for minimal median coverage to keep a zone (float),  [default: %default]')
	parser.add_option( '', '--min_zone', dest='min_zone', default='500', help='Minimal number of covered sites in a zone to be considered (integer),  [default: %default]')
	parser.add_option( '', '--min_gap', dest='min_gap', default='300', help='Maximal number of contiguous uncovered sites in a zone to be considered as a single zone (integer),  [default: %default]')
	parser.add_option( '', '--thread', dest='thread', default='1', help='The thread number used for mapping (integer), [default: %default]')
	parser.add_option( '', '--msd', dest='msd', default='3', help='Multiplicator of standard deviation for discordant zone identification (integer), [default: %default]')
	parser.add_option( '', '--max_dist_merge', dest='max_dist_merge', default=1000, help='Maximal distance between two discordant zone to merge, [default: %default]')
	parser.add_option( '', '--YiS', dest='YiS', default=0, help='The Y-intercept of the linear function for zone size that will give the first component of product giving the score (integer), [default: %default]')
	parser.add_option( '', '--MiS', dest='MiS', default=0.5, help='Multiplicator of median insert size for calculating minimal zone size for which the first component of product giving the score will be maximal (integer), [default: %default]. Exmple: if 0.5, discordant zone of more than 2500 pb will have a maximal score')
	parser.add_option( '', '--YiC', dest='YiC', default=0, help='The Y-intercept of the linear function for coverage that will give the second component of product giving the score (integer), [default: %default]')
	parser.add_option( '', '--MiC', dest='MiC', default=0.25, help='Multiplicator of median coverage for calculating  minimal zone coverage for which the second component of product giving the score will be maximal (integer), [default: %default]. For homozygous SV in diploid: expected value = 0.5, if heterozygous: expected value = 0.25')
	parser.add_option( '', '--min_score', dest='min_score', default=70, help='The minimal score for a discordant zone to be identified as passed, [default: %default]')
	parser.add_option( '', '--ploid', dest='ploid', default=0.33, help='Multiplicator for coverage variation detection in SV identification (ex : If homozygous duplication expected in diploid: expected = coverage + coverage*1, if heterozygous duplication expected in diploid => expected = coverage + coverage*0.5). Choose a value lower than the expected one')
	parser.add_option( '', '--restimate', dest='restimate', default='n', help='Wether re-estimating --mini and --maxi parameters: y or n, [default: %default]. If y, these parameters are calculated as followed on well mapped paired read on the basis of previous min and max parameters: min/max = median -/+ (standard_deviation * "--msd" option)')
	# parser.add_option( '', '--output', dest='output', default='config.conf', help='The output of the conf file, [default: %default]')
	# parser.add_option( '', '--chr', dest='chr', default='chr.tab', help='Output file containing chromosomes informations, [default: %default]')
	parser.add_option( '', '--rm_intermediate', dest='rm_intermediate', default='y', help='Remove intermediate bam/sam, [default: %default]')
	parser.add_option( '', '--prefix', dest='prefix', default='apmap', help='Prefix for all output files, [default: %default]')
	parser.add_option( '', '--step', dest='step', default='1234567', help='Steps to perform, [default: %default]')
	parser.add_option( '', '--exclude_chrom', dest='exclude_chrom', default='no_exclude', help='Exclude chromosomes from analysis. "no_exclude" or chromosomes names separated by "=", [default: %default]')
	(options, args) = parser.parse_args()



	nbProcs = int(options.thread)
	if nbProcs > multiprocessing.cpu_count():
		sys.exit("Processors number to high.\nYou have only "+str(multiprocessing.cpu_count())+" processor(s) available on this computer.")

	pathname = os.path.dirname(os.path.abspath(sys.argv[0]))

	loca_programs = ConfigParser.RawConfigParser()
	loca_programs.read(pathname+'/loca_programs.conf')

	if '1' in options.step:
		t0 = datetime.datetime.now()
		print("Step 1 'create_conf' in progress")
		sys.stdout.flush()
		conf_commande = '%s %s/1_create_conf.py --tool %s --ref %s --q1 %s --q2 %s --orient %s --mini %s --maxi %s --qual %s --index %s --rmindex %s --mini_dis %s --mult_max_cov %s --mult_min_cov %s --min_zone %s --min_gap %s --thread %s --msd %s --max_dist_merge %s --YiS %s --MiS %s --YiC %s --MiC %s --min_score %s --ploid %s --restimate %s --output %s.conf --chr %s.chrom --rm_intermediate %s --exclude_chrom %s' % (loca_programs.get('Programs','python'), pathname, options.tool, options.ref, options.q1, options.q2, options.orient, options.mini, options.maxi, options.qual, options.index, options.rmindex, options.mini_dis, options.mult_max_cov, options.mult_min_cov, options.min_zone, options.min_gap, options.thread, options.msd, options.max_dist_merge, options.YiS, options.MiS, options.YiC, options.MiC, options.min_score, options.ploid, options.restimate, options.prefix, options.prefix, options.rm_intermediate, options.exclude_chrom)
		# print conf_commande
		run_job( conf_commande, 'bug')
		print("Step 1 is finished (time : "+str(datetime.datetime.now()-t0)+")")
		sys.stdout.flush()
	if '2' in options.step:
		t0 = datetime.datetime.now()
		print("Step 2 'map' in progress")
		sys.stdout.flush()
		mapping = '%s %s/2_map.py --config %s.conf --out %s.sam' % (loca_programs.get('Programs','python'), pathname, options.prefix, options.prefix)
		# print mapping
		run_job( mapping, 'bug')
		print("Step 2 is finished (time : "+str(datetime.datetime.now()-t0)+")")
		sys.stdout.flush()
	if '3' in options.step:
		t0 = datetime.datetime.now()
		print("Step 3 'filter_single_pair' in progress")
		sys.stdout.flush()
		filter1 = '%s %s/3_filter_single_pair.py --config %s.conf --out %s_fltr1.sam' % (loca_programs.get('Programs','python'), pathname, options.prefix, options.prefix)
		# print filter1
		run_job( filter1, 'bug')
		print("Step 3 is finished (time : "+str(datetime.datetime.now()-t0)+")")
		sys.stdout.flush()
	if '4' in options.step:
		t0 = datetime.datetime.now()
		print("Step 4 'filter_sam' in progress")
		sys.stdout.flush()
		filter2 = '%s %s/4_filter_sam.py --config %s.conf --out %s_fltr2.bam' % (loca_programs.get('Programs','python'), pathname, options.prefix, options.prefix)
		# print filter2
		run_job( filter2, 'bug')
		print("Step 4 is finished (time : "+str(datetime.datetime.now()-t0)+")")
		sys.stdout.flush()
	if '5' in options.step:
		t0 = datetime.datetime.now()
		print("Step 5 'calc_stat' in progress")
		sys.stdout.flush()
		stat = '%s %s/5_calc_stat.py --config %s.conf --out %s.cov --stat %s_stat.txt --outconf %s_stat.conf' % (loca_programs.get('Programs','python'), pathname, options.prefix, options.prefix, options.prefix, options.prefix)
		# print stat
		run_job( stat, 'bug')
		print("Step 5 is finished (time : "+str(datetime.datetime.now()-t0)+")")
		sys.stdout.flush()
	if '6' in options.step:
		t0 = datetime.datetime.now()
		print("Step 6 'parse_discord' in progress")
		sys.stdout.flush()
		trie = '%s %s/6_parse_discord.py --config %s.conf --out_ins %s_ins.bam --out_del %s_del.bam --out_fr %s_fr.bam --out_rf %s_rf.bam --out_ff %s_ff.bam --out_rr %s_rr.bam --out_chr_fr %s_chr_fr.bam --out_chr_rf %s_chr_rf.bam --out_chr_ff %s_chr_ff.bam --out_chr_rr %s_chr_rr.bam --out_discarded %s_discarded.bam --liste_type %s.list --discord_prop %s.prop' % (loca_programs.get('Programs','python'), pathname, options.prefix, options.prefix, options.prefix, options.prefix, options.prefix, options.prefix, options.prefix, options.prefix, options.prefix, options.prefix, options.prefix, options.prefix, options.prefix, options.prefix)
		# print trie
		run_job( trie, 'bug')
		print("Step 6 is finished (time : "+str(datetime.datetime.now()-t0)+")")
		sys.stdout.flush()
	if '7' in options.step:
		t0 = datetime.datetime.now()
		print("Step 7 'select_on_cov' in progress")
		sys.stdout.flush()
		liste_id = []
		if options.orient == 'rf':
			select1 = '%s %s/7_select_on_cov.py --config %s.conf --sam %s_chr_fr.bam --out %s_chr_fr.score' % (loca_programs.get('Programs','python'), pathname, options.prefix, options.prefix, options.prefix)
			liste_id.append(select1)
			select2 = '%s %s/7_select_on_cov.py --config %s.conf --sam %s_chr_ff.bam --out %s_chr_ff.score' % (loca_programs.get('Programs','python'), pathname, options.prefix, options.prefix, options.prefix)
			liste_id.append(select2)
			select3 = '%s %s/7_select_on_cov.py --config %s.conf --sam %s_chr_rr.bam --out %s_chr_rr.score' % (loca_programs.get('Programs','python'), pathname, options.prefix, options.prefix, options.prefix)
			liste_id.append(select3)
			select4 = '%s %s/7_select_on_cov.py --config %s.conf --sam %s_chr_rf.bam --out %s_chr_rf.score' % (loca_programs.get('Programs','python'), pathname, options.prefix, options.prefix, options.prefix)
			liste_id.append(select4)
			select5 = '%s %s/7_select_on_cov.py --config %s.conf --sam %s_del.bam --out %s_del.score' % (loca_programs.get('Programs','python'), pathname, options.prefix, options.prefix, options.prefix)
			liste_id.append(select5)
			select6 = '%s %s/7_select_on_cov.py --config %s.conf --sam %s_ins.bam --out %s_ins.score' % (loca_programs.get('Programs','python'), pathname, options.prefix, options.prefix, options.prefix)
			liste_id.append(select6)
			select7 = '%s %s/7_select_on_cov.py --config %s.conf --sam %s_fr.bam --out %s_fr.score' % (loca_programs.get('Programs','python'), pathname, options.prefix, options.prefix, options.prefix)
			liste_id.append(select7)
			select8 = '%s %s/7_select_on_cov.py --config %s.conf --sam %s_ff.bam --out %s_ff.score' % (loca_programs.get('Programs','python'), pathname, options.prefix, options.prefix, options.prefix)
			liste_id.append(select8)
			select9 = '%s %s/7_select_on_cov.py --config %s.conf --sam %s_rr.bam --out %s_rr.score' % (loca_programs.get('Programs','python'), pathname, options.prefix, options.prefix, options.prefix)
			liste_id.append(select9)


			#Change the chromosome file to working in different sub folders
			config = ConfigParser.RawConfigParser()
			config.read(options.prefix+".conf")
			chrFile = config.get('General','chr')
			config.set('General','chr', "../"+chrFile)
			with open(options.prefix+".conf", 'wb') as configfile:
				config.write(configfile)

			pool = multiprocessing.Pool(processes=nbProcs)
			resultsJobs = pool.map(main, liste_id)

			for i, job in enumerate(resultsJobs):
				if job != 0:
					print("Sorry the job : \n"+liste_id[i]+"\n" \
						  "could not be completed due to an error.\n" \
						  "Please read the error log file for more details.")

			#rewrite the chromosome file to his initial value
			chrFile = config.get('General','chr')[3:]
			config.set('General','chr', chrFile)
			with open(options.prefix+".conf", 'wb') as configfile:
				config.write(configfile)

		elif options.orient == 'fr':
			select1 = '%s %s/7_select_on_cov.py --config %s.conf --sam %s_chr_fr.bam --out %s_chr_fr.score' % (loca_programs.get('Programs','python'), pathname, options.prefix, options.prefix, options.prefix)
			liste_id.append(select1)
			select2 = '%s %s/7_select_on_cov.py --config %s.conf --sam %s_chr_ff.bam --out %s_chr_ff.score' % (loca_programs.get('Programs','python'), pathname, options.prefix, options.prefix, options.prefix)
			liste_id.append(select2)
			select3 = '%s %s/7_select_on_cov.py --config %s.conf --sam %s_chr_rr.bam --out %s_chr_rr.score' % (loca_programs.get('Programs','python'), pathname, options.prefix, options.prefix, options.prefix)
			liste_id.append(select3)
			select4 = '%s %s/7_select_on_cov.py --config %s.conf --sam %s_chr_rf.bam --out %s_chr_rf.score' % (loca_programs.get('Programs','python'), pathname, options.prefix, options.prefix, options.prefix)
			liste_id.append(select4)
			select5 = '%s %s/7_select_on_cov.py --config %s.conf --sam %s_del.bam --out %s_del.score' % (loca_programs.get('Programs','python'), pathname, options.prefix, options.prefix, options.prefix)
			liste_id.append(select5)
			select6 = '%s %s/7_select_on_cov.py --config %s.conf --sam %s_ins.bam --out %s_ins.score' % (loca_programs.get('Programs','python'), pathname, options.prefix, options.prefix, options.prefix)
			liste_id.append(select6)
			select7 = '%s %s/7_select_on_cov.py --config %s.conf --sam %s_rf.bam --out %s_rf.score' % (loca_programs.get('Programs','python'), pathname, options.prefix, options.prefix, options.prefix)
			liste_id.append(select7)
			select8 = '%s %s/7_select_on_cov.py --config %s.conf --sam %s_ff.bam --out %s_ff.score' % (loca_programs.get('Programs','python'), pathname, options.prefix, options.prefix, options.prefix)
			liste_id.append(select8)
			select9 = '%s %s/7_select_on_cov.py --config %s.conf --sam %s_rr.bam --out %s_rr.score' % (loca_programs.get('Programs','python'), pathname, options.prefix, options.prefix, options.prefix)
			liste_id.append(select9)


			#Change the chromosome file to working in different sub folders
			config = ConfigParser.RawConfigParser()
			config.read(options.prefix+".conf")
			chrFile = config.get('General','chr')
			config.set('General','chr', "../"+chrFile)
			with open(options.prefix+".conf", 'wb') as configfile:
				config.write(configfile)

			pool = multiprocessing.Pool(processes=nbProcs)
			resultsJobs = pool.map(main, liste_id)

			for i, job in enumerate(resultsJobs):
				if job != 0:
					print("Sorry the job : \n"+liste_id[i]+"\n" \
						  "could not be completed due to an error.\n" \
						  "Please read the error log file for more details.")

			#rewrite the chromosome file to his initial value
			chrFile = config.get('General','chr')[3:]
			config.set('General','chr', chrFile)
			with open(options.prefix+".conf", 'wb') as configfile:
				config.write(configfile)
		else:
			mot = 'Unrecognized argument in --orient option: '+options.orient
			sys.exit(mot)
		print("Step 7 is finished (time : "+str(datetime.datetime.now()-t0)+")")
		sys.stdout.flush()
if __name__ == "__main__": __main__()
