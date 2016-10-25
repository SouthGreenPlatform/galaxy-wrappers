
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



def __main__():
	#Parse Command Line
	parser = optparse.OptionParser(usage="python %prog [options]\n\nProgram designed by Guillaume MARTIN : guillaume.martin@cirad.fr\n\n"
	"This script takes scaffold to join and group them by linkage. The input is a tabulated file looking as folowed:\n"
	"scaffold1    1       2458028 FWD     scaffold2    1       1       FWD     FWD     contig\n"
	"scaffold36   1       250001  FWD     scaffold5    2000000 2000000 FWD     REV     contig")
	# Wrapper options. 
	parser.add_option( '', '--table', dest='table', help='The table file input')
	parser.add_option( '', '--out', dest='out', default='intermediate_junction.txt', help='The output file name, [default: %default]')
	(options, args) = parser.parse_args()
	
	# Filling dictionnary
	dico_line = set()
	file = open(options.table)
	for line in file:
		data = line.split()
		if data:
			if data[3] != 'FWD' or data[7] != 'FWD':
				sys.exit('Warning! the orientation filled in column 4 or 8 is not managed, this column should contain "FWD" ')
			dico_line.add(line)
	file.close()
	
	dico = set() # --> to identify already treated lines
	dico_scaff = {} # --> a hash table containing a list of three elements : (1) grouped scaffold ordered, (2) scaffold orientation in the group, (3) ungrouped scaffolds, (4) scaffold number
	i = 0
	for line in dico_line:
		data = line.split()
		if data:
			if not(line in dico): # --> this line start a new group
				i += 1
				# print i
				dico.add(line) # --> record a line already treated
				dico_scaff[i] = [[data[4]],[data[7]],set(), 1]
				taille_dico = int(dico_scaff[i][3]) # --> to define a group is inflating or not (in the while loop)
				if data[5] == '1':# --> scaffold should be inserted at the begining of the list
					dico_scaff[i][0].insert(0,data[0])
					dico_scaff[i][1].insert(0,data[8])
					dico_scaff[i][3] += 1
				else:
					dico_scaff[i][0].append(data[0])
					dico_scaff[i][1].append(data[8])
					dico_scaff[i][3] += 1
				# --> The new group as been stared
				while taille_dico < dico_scaff[i][3]: # --> verification is new scaffolds have been added to group
					taille_dico = int(dico_scaff[i][3])
					# print dico_scaff[i][0], dico_scaff[i][1], len(dico_scaff[i][2]), dico_scaff[i][3]
					for line2 in dico_line:
						data2 = line2.split()
						if data2 and not(line2 in dico):
							if data2[4] in dico_scaff[i][0]: # --> checking if this scaffold is already in the group
								dico.add(line2)
								dico_scaff[i][3] += 1
								if data2[4] == dico_scaff[i][0][-1]: # --> this scaffold is at the end of the list
									if dico_scaff[i][1][-1] == 'FWD': # --> if this scaffold is FWD orientated
										if data2[5] == '1': # --> the new scaffold could not be inserted
											dico_scaff[i][2].add(line2)
											# print 'toto1'
										else: # --> the new scaffold could be inserted at the end
											dico_scaff[i][0].append(data2[0])
											dico_scaff[i][1].append(data2[8])
									elif dico_scaff[i][1][-1] == 'REV': # --> if this scaffold is REV orientated
										if data2[5] != '1': # --> the new scaffold could not be inserted
											dico_scaff[i][2].add(line2)
											# print 'toto2'
										else: # --> the new scaffold could be inserted at the end but its orientation should be reverted
											dico_scaff[i][0].append(data2[0])
											if data2[8] == 'FWD':
												dico_scaff[i][1].append('REV')
											elif data2[8] == 'REV':
												dico_scaff[i][1].append('FWD')
											else:
												sys.exit('Unrecognized orientaion in column 8. Only FWD and REV are recognized')
									else:
										sys.exit('Unrecognized orientaion in column 8. Only FWD and REV are recognized')
								elif data2[4] == dico_scaff[i][0][0]: # --> this scaffold is at the begining of the list
									if dico_scaff[i][1][0] == 'FWD': # --> if this scaffold is FWD orientated
										if data2[5] == '1':  # --> the new scaffold could be inserted at the begining
											dico_scaff[i][0].insert(0,data2[0])
											dico_scaff[i][1].insert(0,data2[8])
										else:# --> the new scaffold could not be inserted
											dico_scaff[i][2].add(line2)
											# print 'toto3'
									elif dico_scaff[i][1][0] == 'REV': # --> if this scaffold is REV orientated
										if data2[5] == '1': # --> the new scaffold could not be inserted
											dico_scaff[i][2].add(line2)
											# print 'toto4'
										else: # --> the new scaffold could be inserted at the begining but its orientation should be reverted
											dico_scaff[i][0].insert(0,data2[0])
											if data2[8] == 'FWD':
												dico_scaff[i][1].insert(0,'REV')
											elif data2[8] == 'REV':
												dico_scaff[i][1].insert(0,'FWD')
											else:
												sys.exit('Unrecognized orientaion in column 8. Only FWD and REV are recognized')
									else:
										sys.exit('Unrecognized orientaion in column 8. Only FWD and REV are recognized')
								else:
									dico_scaff[i][2].add(line2)
									# print 'toto5'
							elif data2[0] in dico_scaff[i][0]: # --> checking if this scaffold is already in the group
								dico.add(line2)
								dico_scaff[i][3] += 1
								if data2[0] == dico_scaff[i][0][-1]: # --> this scaffold is at the end of the list
									if dico_scaff[i][1][-1] == 'FWD': # --> if this scaffold is FWD orientated
										if data2[5] == '1':
											if data2[8] == 'FWD': # --> the new scaffold could be inserted at the end
												dico_scaff[i][0].append(data2[4])
												dico_scaff[i][1].append(data2[7])
											elif data2[8] == 'REV': # --> the new scaffold could not be inserted
												dico_scaff[i][2].add(line2)
												# print 'toto6'
											else:
												sys.exit('Unrecognized orientaion in column 8. Only FWD and REV are recognized')
										else:
											if data2[8] == 'REV': # --> the new scaffold could be inserted at the end
												dico_scaff[i][0].append(data2[4])
												dico_scaff[i][1].append('REV')
											elif data2[8] == 'FWD': # --> the new scaffold could not be inserted
												dico_scaff[i][2].add(line2)
												# print 'toto7'
											else:
												sys.exit('Unrecognized orientaion in column 8. Only FWD and REV are recognized')
									elif dico_scaff[i][1][-1] == 'REV': # --> if this scaffold is REV orientated
										if data2[5] == '1':
											if data2[8] == 'REV': # --> the new scaffold could be inserted at the end
												dico_scaff[i][0].append(data2[4])
												dico_scaff[i][1].append(data2[7])
											elif data2[8] == 'FWD': # --> the new scaffold could not be inserted
												dico_scaff[i][2].add(line2)
												# print 'toto8'
											else:
												sys.exit('Unrecognized orientaion in column 8. Only FWD and REV are recognized')
										else:
											if data2[8] == 'FWD': # --> the new scaffold could be inserted at the end
												dico_scaff[i][0].append(data2[4])
												dico_scaff[i][1].append('REV')
											elif data2[8] == 'REV': # --> the new scaffold could not be inserted
												dico_scaff[i][2].add(line2)
												# print 'toto9'
											else:
												sys.exit('Unrecognized orientaion in column 8. Only FWD and REV are recognized')
									else:
										sys.exit('Unrecognized orientaion in column 8. Only FWD and REV are recognized')
								elif data2[0] == dico_scaff[i][0][0]: # --> this scaffold is at the begining of the list
									if dico_scaff[i][1][0] == 'FWD': # --> if this scaffold is FWD orientated
										if data2[5] != '1':
											if data2[8] == 'FWD': # --> the new scaffold could be inserted at the begining
												dico_scaff[i][0].insert(0,data2[4])
												dico_scaff[i][1].insert(0,data2[7])
											elif data2[8] == 'REV': # --> the new scaffold could not be inserted
												dico_scaff[i][2].add(line2)
												# print 'toto10'
											else:
												sys.exit('Unrecognized orientaion in column 8. Only FWD and REV are recognized')
										else:
											if data2[8] == 'REV': # --> the new scaffold could be inserted at begining
												dico_scaff[i][0].insert(0,data2[4])
												dico_scaff[i][1].insert(0,'REV')
											elif data2[8] == 'FWD': # --> the new scaffold could not be inserted
												dico_scaff[i][2].add(line2)
												# print 'toto11'
											else:
												sys.exit('Unrecognized orientaion in column 8. Only FWD and REV are recognized')
									elif dico_scaff[i][1][0] == 'REV': # --> if this scaffold is REV orientated
										if data2[5] == '1':
											if data2[8] == 'FWD': # --> the new scaffold could be inserted at the begining
												dico_scaff[i][0].insert(0,data2[4])
												dico_scaff[i][1].insert(0,'REV')
											elif data2[8] == 'REV': # --> the new scaffold could not be inserted
												dico_scaff[i][2].add(line2)
												# print 'toto12'
											else:
												sys.exit('Unrecognized orientaion in column 8. Only FWD and REV are recognized')
										else:
											if data2[8] == 'REV': # --> the new scaffold could be inserted at begining
												dico_scaff[i][0].insert(0,data2[4])
												dico_scaff[i][1].insert(0,data2[7])
											elif data2[8] == 'FWD': # --> the new scaffold could not be inserted
												dico_scaff[i][2].add(line2)
												# print 'toto13'
											else:
												sys.exit('Unrecognized orientaion in column 8. Only FWD and REV are recognized')
									else:
										sys.exit('Unrecognized orientaion in column 8. Only FWD and REV are recognized')
								else:
									dico_scaff[i][2].add(line2)
									# print 'toto14'

	
	outfile = open(options.out,'w')
	for n in dico_scaff:
		liste = list(dico_scaff[n][0])
		liste.sort()
		outfile.write('>'+liste[0]+'\n')
		# print '>'+liste[0]
		for k in dico_scaff[n][0]:
			outfile.write('\t'.join([k,dico_scaff[n][1][dico_scaff[n][0].index(k)]])+'\n')
			# print '\t'.join([k,dico_scaff[n][1][dico_scaff[n][0].index(k)]])
		for k in dico_scaff[n][2]:
			data = k.split()
			if data[5] == '1':
				if data[3] == data[8]:
					outfile.write('not_grouped\t'+data[0]+'\tFWD\t'+data[4]+'\t'+data[7]+'\n')
				else:
					outfile.write('not_grouped\t'+data[0]+'\tREV\t'+data[4]+'\t'+data[7]+'\n')
			else:
				if data[3] == data[8]:
					outfile.write('not_grouped\t'+data[4]+'\t'+data[7]+'\t'+data[0]+'\tFWD'+'\n')
				else:
					outfile.write('not_grouped\t'+data[4]+'\t'+data[7]+'\t'+data[0]+'\tREV'+'\n')
	outfile.close()
	
	
	
if __name__ == "__main__": __main__()