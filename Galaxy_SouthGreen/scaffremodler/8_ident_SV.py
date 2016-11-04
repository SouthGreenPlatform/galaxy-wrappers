#!/usr/local/bioinfo/python/2.7.9/bin/python
#
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

import optparse
import os
import shutil
import subprocess
import sys
import tempfile
import fileinput
import ConfigParser
import operator
import time
import random
import datetime
import ctypes
import multiprocessing as mp
from multiprocessing.sharedctypes import Value, Array 

# Global variables
covChr = {}  # contains the coverture of each site by chromosome


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


def readChrLenght(chrFile):
	"""
		Return a dict with the chromosome name as key, and their size in value.
	"""
	chrSize = {}
	f = open(chrFile, 'r')
	for line in f:
		if line.strip():
			cols = line.split()
			chrSize[cols[0]] = int(cols[1])
	f.close()
	return chrSize

	
def readCov(chrFile, covFile):
	"""
		fill a global dict with the coverage of the chromosome for each site
	"""
	chrSize = readChrLenght(chrFile)
	chr = ""
	f = open(covFile, 'r')
	for line in f:
		if line.strip():
			cols = line.split()
			if cols[0] != chr:  # new chromosome
				if chr:  # not the first line
					covChr[chr] = test
				test = mp.sharedctypes.RawArray(ctypes.c_int, chrSize[cols[0]])
				test[int(cols[1])-1] = int(cols[2])
				chr = cols[0]
			else:
				test[int(cols[1])-1] = int(cols[2])
	covChr[chr] = test

	
def indent_discord(FF, FR, RR, INS, DEL, CHR_rr, CHR_fr, CHR_rf, CHR_ff, INSERT, OUT, EXP_COV, PLOID, TYPE):
	outfile = open(OUT,'w')
	file = open(FF)
	dic_FF = {}
	for line in file:
		data = line.split()
		if data != []:
			if data[0][0] != "#" and data[13] == "PASSED":
				if data[0] in dic_FF:
					dic_FF[data[0]].append([int(data[1]), int(data[2]), data[5], int(data[6]), int(data[7])])
				else :
					dic_FF[data[0]] = ([])
					dic_FF[data[0]].append([int(data[1]), int(data[2]), data[5], int(data[6]), int(data[7])])
	file.close()
	
	file = open(FR)
	dic_FR = {}
	for line in file:
		data = line.split()
		if data != []:
			if data[0][0] != "#" and data[13] == "PASSED":
				if data[0] in dic_FR:
					dic_FR[data[0]].append([int(data[1]), int(data[2]), data[5], int(data[6]), int(data[7])])
				else :
					dic_FR[data[0]] = ([])
					dic_FR[data[0]].append([int(data[1]), int(data[2]), data[5], int(data[6]), int(data[7])])
	file.close()
	
	file = open(RR)
	dic_RR = {}
	for line in file:
		data = line.split()
		if data != []:
			if data[0][0] != "#" and data[13] == "PASSED":
				if data[0] in dic_RR:
					dic_RR[data[0]].append([int(data[1]), int(data[2]), data[5], int(data[6]), int(data[7])])
				else :
					dic_RR[data[0]] = ([])
					dic_RR[data[0]].append([int(data[1]), int(data[2]), data[5], int(data[6]), int(data[7])])
	file.close()
	
	file = open(INS)
	dic_INS = {}
	for line in file:
		data = line.split()
		if data != []:
			if data[0][0] != "#" and data[13] == "PASSED":
				if data[0] in dic_INS:
					dic_INS[data[0]].append([int(data[1]), int(data[2]), data[5], int(data[6]), int(data[7])])
				else :
					dic_INS[data[0]] = ([])
					dic_INS[data[0]].append([int(data[1]), int(data[2]), data[5], int(data[6]), int(data[7])])
	file.close()
	
	file = open(DEL)
	dic_DEL = {}
	for line in file:
		data = line.split()
		if data != []:
			if data[0][0] != "#" and data[13] == "PASSED":
				if data[0] in dic_DEL:
					dic_DEL[data[0]].append([int(data[1]), int(data[2]), data[5], int(data[6]), int(data[7])])
				else :
					dic_DEL[data[0]] = ([])
					dic_DEL[data[0]].append([int(data[1]), int(data[2]), data[5], int(data[6]), int(data[7])])
	file.close()
	
	file = open(CHR_rr)
	dic_CHR_rr = {}
	dic_CHR_rr_rev = {}
	for line in file:
		data = line.split()
		if data != []:
			if data[0][0] != "#" and data[13] == "PASSED":
				if data[0] in dic_CHR_rr:
					dic_CHR_rr[data[0]].append([int(data[1]), int(data[2]), data[5], int(data[6]), int(data[7])])
				else :
					dic_CHR_rr[data[0]] = ([])
					dic_CHR_rr[data[0]].append([int(data[1]), int(data[2]), data[5], int(data[6]), int(data[7])])
				if data[5] in dic_CHR_rr_rev:
					dic_CHR_rr_rev[data[5]].append([int(data[6]), int(data[7]), data[0], int(data[1]), int(data[2])])
				else :
					dic_CHR_rr_rev[data[5]] = ([])
					dic_CHR_rr_rev[data[5]].append([int(data[6]), int(data[7]), data[0], int(data[1]), int(data[2])])
	file.close()
	
	
	file = open(CHR_rf)
	dic_CHR_rf = {}
	dic_CHR_rf_rev = {}
	for line in file:
		data = line.split()
		if data != []:
			if data[0][0] != "#" and data[13] == "PASSED":
				if data[0] in dic_CHR_rf:
					dic_CHR_rf[data[0]].append([int(data[1]), int(data[2]), data[5], int(data[6]), int(data[7])])
				else :
					dic_CHR_rf[data[0]] = ([])
					dic_CHR_rf[data[0]].append([int(data[1]), int(data[2]), data[5], int(data[6]), int(data[7])])
				if data[5] in dic_CHR_rf_rev:
					dic_CHR_rf_rev[data[5]].append([int(data[6]), int(data[7]), data[0], int(data[1]), int(data[2])])
				else :
					dic_CHR_rf_rev[data[5]] = ([])
					dic_CHR_rf_rev[data[5]].append([int(data[6]), int(data[7]), data[0], int(data[1]), int(data[2])])
	file.close()
	
	file = open(CHR_fr)
	dic_CHR_fr = {}
	dic_CHR_fr_rev = {}
	for line in file:
		data = line.split()
		if data != []:
			if data[0][0] != "#" and data[13] == "PASSED":
				if data[0] in dic_CHR_fr:
					dic_CHR_fr[data[0]].append([int(data[1]), int(data[2]), data[5], int(data[6]), int(data[7])])
				else :
					dic_CHR_fr[data[0]] = ([])
					dic_CHR_fr[data[0]].append([int(data[1]), int(data[2]), data[5], int(data[6]), int(data[7])])
				if data[5] in dic_CHR_fr_rev:
					dic_CHR_fr_rev[data[5]].append([int(data[6]), int(data[7]), data[0], int(data[1]), int(data[2])])
				else :
					dic_CHR_fr_rev[data[5]] = ([])
					dic_CHR_fr_rev[data[5]].append([int(data[6]), int(data[7]), data[0], int(data[1]), int(data[2])])
	file.close()
	
	file = open(CHR_ff)
	dic_CHR_ff = {}
	dic_CHR_ff_rev = {}
	for line in file:
		data = line.split()
		if data != []:
			if data[0][0] != "#" and data[13] == "PASSED":
				if data[0] in dic_CHR_ff:
					dic_CHR_ff[data[0]].append([int(data[1]), int(data[2]), data[5], int(data[6]), int(data[7])])
				else :
					dic_CHR_ff[data[0]] = ([])
					dic_CHR_ff[data[0]].append([int(data[1]), int(data[2]), data[5], int(data[6]), int(data[7])])
				if data[5] in dic_CHR_ff_rev:
					dic_CHR_ff_rev[data[5]].append([int(data[6]), int(data[7]), data[0], int(data[1]), int(data[2])])
				else :
					dic_CHR_ff_rev[data[5]] = ([])
					dic_CHR_ff_rev[data[5]].append([int(data[6]), int(data[7]), data[0], int(data[1]), int(data[2])])
	file.close()
	vasouille = 100
	#All discordant position has been loaded
	##################################################################
	#####This part only works for intra chromosomal rearrangments#####
	##################################################################
	#Now it's time to detect simple translocations
	if TYPE == '0':
		for n in dic_FR:
			for j in dic_FR[n]:
				pos1 = [j[0], j[1]]
				pos3 = [j[3], j[4]]
				pos2 = []
				pos4 = []
				if n in dic_DEL:
					for k in dic_DEL[n]:
						if k[3] <= pos3[1] + 2*INSERT and  k[3] + vasouille >= pos3[1] and k[0] + vasouille >= pos1[1] and k[1] <= pos3[0] + vasouille:#potential translocation
							pos2 = [k[0], k[1]]
							pos4 = [k[3], k[4]]
							for m in dic_DEL[n]:
								if m[1] <= pos1[0] + vasouille and m[1] + 2*INSERT >= pos1[0] and m[3] + vasouille >= pos2[1] and m[3] <= pos2[1] + 2*INSERT:
									pos0 = [m[0],m[1]]
									pos5 = [m[3],m[4]]
									outfile.write('translocation\tregion:\t'+n+'\t'+str(pos5[0])+'\t'+str(pos3[1])+'\ttarget:\t'+n+'\t'+str(pos0[1])+'\t'+str(pos1[0])+'\tALTERNATIVE\tregion:\t'+n+'\t'+str(pos1[0])+'\t'+str(pos2[1])+'\ttarget:\t'+n+'\t'+str(pos3[1])+'\t'+str(pos4[0])+'\n')
	#Now it's time to detect simple duplications
	if TYPE == '1':
		for n in dic_FR:
			for j in dic_FR[n]:
				pos1 = [j[0], j[1]]
				pos3 = [j[3], j[4]]
				pos2 = []
				pos4 = []
				if n in dic_DEL:
					for k in dic_DEL[n]:
						if k[3] <= pos3[1] + 2*INSERT and  k[3] + vasouille >= pos3[1] and k[0] + vasouille >= pos1[1] and k[1] <= pos3[0] + vasouille:#potential duplication
							pos2 = [k[0], k[1]]
							pos4 = [k[3], k[4]]
							#Maintenant s'occuper des couvertures
							if couv_med_reg(pos1[0], pos2[1], n) >= (EXP_COV + float(EXP_COV)*PLOID):
								outfile.write('duplication\tregion:\t'+n+'\t'+str(pos1[0])+'\t'+str(pos2[1])+'\ttarget:\t'+n+'\t'+str(pos3[1])+'\t'+str(pos4[0])+'\n')
						elif pos1[0] <= k[1] + 2*INSERT and  pos1[0] + vasouille >= k[1] and pos3[0] + vasouille >= k[4] and pos1[1] <= k[3] + vasouille:#potential duplication
							pos2 = [k[0], k[1]]
							pos4 = [k[3], k[4]]
							#Maintenant s'occuper des couvertures
							if couv_med_reg(pos4[0], pos3[1], n) >= (EXP_COV + float(EXP_COV)*PLOID):
								outfile.write('duplication\tregion:\t'+n+'\t'+str(pos4[0])+'\t'+str(pos3[1])+'\ttarget:\t'+n+'\t'+str(pos2[1])+'\t'+str(pos1[0])+'\n')
	#Now it's time to detect simple inversions
	if TYPE == '0':
		for n in dic_RR:
			for j in dic_RR[n]:
				pos1 = [j[0], j[1]]
				pos3 = [j[3], j[4]]
				pos2 = []
				pos4 = []
				if n in dic_FF:
					for k in dic_FF[n]:
						if k[0] <= pos1[1] + 2*INSERT and  k[0] + vasouille >= pos1[1] and k[3] <= pos3[1] + 2*INSERT and k[3] + vasouille >= pos3[1]:#simple inversion
							pos2 = [k[0], k[1]]
							pos4 = [k[3], k[4]]
							outfile.write('inversion\tregion:\t'+n+'\t'+str(pos2[0])+'\t'+str(pos3[1])+'\n')
	#Now it's time to detect simple deletions
	if TYPE == '2':
		for n in dic_DEL:
			for j in dic_DEL[n]:
				pos1 = [j[0], j[1]]
				pos3 = [j[3], j[4]]
				if couv_med_reg_del(pos1[1], pos3[0], n) <= (EXP_COV - float(EXP_COV)*PLOID):
					outfile.write('deletion\tregion:\t'+n+'\t'+str(pos1[1])+'\t'+str(pos3[0])+'\n')
	#Now it's time to detect simple tandem duplication
	if TYPE == '3':
		for n in dic_FR:
			for j in dic_FR[n]:
				pos1 = [j[0], j[1]]
				pos3 = [j[3], j[4]]
				if couv_med_reg(pos1[0], pos3[1], n) >= (EXP_COV + float(EXP_COV)*PLOID):
					outfile.write('tandem_duplication\tregion:\t'+n+'\t'+str(pos1[0])+'\t'+str(pos3[1])+'\n')
	#Now it's time to detect simple reciprocal translocations
	if TYPE == '0':
		for n in dic_FR:
			dic_DEL_FR = []
			for j in dic_FR[n]:
				pos1 = [j[0], j[1]]
				pos3 = [j[3], j[4]]
				if n in dic_DEL:
					for k in dic_DEL[n]:
						if pos1[0] <= k[1] + 2*INSERT and pos1[0] + vasouille >= k[1] and k[3] <= pos3[1] + 2*INSERT and k[3] + vasouille >= pos3[1]:
							pos2 = [k[0], k[1]]
							pos4 = [k[3], k[4]]
							dic_DEL_FR.append([pos2[0], pos1[1], pos3[0], pos4[1], pos1[0], pos4[0], pos2[1], pos3[1]])
			L_done = []
			for x in dic_DEL_FR:
				L_done.append(x)
				for y in dic_DEL_FR:
					if not(y in L_done):
						if x[1] <= y[0] + vasouille and y[1] <= x[2] + vasouille and x[3] <= y[2] + vasouille:
							outfile.write('reciprocal_translocation\tregion1:\t'+n+'\t'+str(x[4])+'\t'+str(y[6])+'\tregion2:\t'+n+'\t'+str(x[5])+'\t'+str(y[7])+'\n')
						elif y[1] <= x[0] + vasouille and x[1] <= y[2] + vasouille and y[3] <= x[2] + vasouille:
							outfile.write('reciprocal_translocation\tregion1:\t'+n+'\t'+str(y[4])+'\t'+str(x[6])+'\tregion2:\t'+n+'\t'+str(y[5])+'\t'+str(x[7])+'\n')
	#Now it's time to detect inversion of a translocation
	if TYPE == '0':
		for n in dic_FF:
			for j in dic_FF[n]:
				pos1 = [j[0], j[1]]
				pos3 = [j[3], j[4]]
				pos2 = []
				pos4 = []
				if n in dic_RR:
					for k in dic_RR[n]:
						if pos3[0] <= k[4] + 2*INSERT and  pos3[0] + vasouille >= k[4] and k[0] + vasouille >= pos1[1]:# F RR F structure
							pos2 = [k[0], k[1]]
							pos4 = [k[3], k[4]]
							if n in dic_DEL:
								for m in dic_DEL[n]:
									if m[1] <= pos1[0] + vasouille and m[1] + 2*INSERT >= pos1[0] and m[3] + vasouille >= pos2[1] and m[3] <= pos2[1] + 2*INSERT:# DEL F R DEL R F structure
										pos0 = [m[0],m[1]]
										pos5 = [m[3],m[4]]
										outfile.write('translocation\tregion_inv:\t'+n+'\t'+str(pos1[0])+'\t'+str(pos2[1])+'\ttarget:\t'+n+'\t'+str(pos4[1])+'\t'+str(pos3[0])+'\n')
		for n in dic_FF:
			for j in dic_FF[n]:
				pos1 = [j[0], j[1]]
				pos3 = [j[3], j[4]]
				pos2 = []
				pos4 = []
				if n in dic_RR:
					for k in dic_RR[n]:
						if pos1[0] <= k[1] + 2*INSERT and  pos1[0] + vasouille >= k[1] and k[3] + vasouille >= pos3[1]:# R FF R structure
							pos2 = [k[0], k[1]]
							pos4 = [k[3], k[4]]
							if n in dic_DEL:
								for m in dic_DEL[n]:
									if m[1] <= pos3[0] + vasouille and m[1] + 2*INSERT >= pos3[0] and m[3] + vasouille >= pos4[1] and m[3] <= pos4[1] + 2*INSERT:# R F DEL F R DEL structure
										pos0 = [m[0],m[1]]
										pos5 = [m[3],m[4]]
										outfile.write('translocation\tregion_inv:\t'+n+'\t'+str(pos3[0])+'\t'+str(pos4[1])+'\ttarget:\t'+n+'\t'+str(pos2[1])+'\t'+str(pos1[0])+'\n')
	#Now it's time to detect inversion of a duplication
	if TYPE == '4':
		for n in dic_FF:
			for j in dic_FF[n]:
				pos1 = [j[0], j[1]]
				pos3 = [j[3], j[4]]
				pos2 = []
				pos4 = []
				if n in dic_RR:
					for k in dic_RR[n]:
						if pos3[0] <= k[4] + 2*INSERT and  pos3[0] + vasouille >= k[4] and k[0] + vasouille >= pos1[1]:# F RR F structure
							pos2 = [k[0], k[1]]
							pos4 = [k[3], k[4]]
							if couv_med_reg(pos1[0], pos2[1], n) >= (EXP_COV + float(EXP_COV)*PLOID):
								outfile.write('duplication\tregion_inv:\t'+n+'\t'+str(pos1[0])+'\t'+str(pos2[1])+'\ttarget:\t'+n+'\t'+str(pos4[1])+'\t'+str(pos3[0])+'\n')
	if TYPE == '5':
		for n in dic_FF:
			for j in dic_FF[n]:
				pos1 = [j[0], j[1]]
				pos3 = [j[3], j[4]]
				pos2 = []
				pos4 = []
				if n in dic_RR:
					for k in dic_RR[n]:
						if pos1[0] <= k[1] + 2*INSERT and  pos1[0] + vasouille >= k[1] and k[3] + vasouille >= pos3[1]:# R FF R structure
							pos2 = [k[0], k[1]]
							pos4 = [k[3], k[4]]
							if couv_med_reg(pos3[0], pos4[1], n) >= (EXP_COV + float(EXP_COV)*PLOID):
								outfile.write('duplication\tregion_inv:\t'+n+'\t'+str(pos3[0])+'\t'+str(pos4[1])+'\ttarget:\t'+n+'\t'+str(pos2[1])+'\t'+str(pos1[0])+'\n')
	#Now it's time to detect inversion of both fragments of reciprocal translocation
	if TYPE == '0':
		for n in dic_FF:
			for j in dic_FF[n]:
				pos1 = [j[0], j[1]]
				pos3 = [j[3], j[4]]
				pos2 = []
				pos4 = []
				if n in dic_RR:
					for k in dic_RR[n]:
						if pos1[0] <= k[1] + 2*INSERT and  pos1[0] + vasouille >= k[1] and pos3[0] <= k[4] + 2*INSERT and  pos3[0] + vasouille >= k[4]:
							pos2 = [k[0], k[1]]
							pos4 = [k[3], k[4]]
							for x in dic_FF[n]:
								if x != j:
									pos5 = [x[0], x[1]]
									pos6 = [x[3], x[4]]
									for y in dic_RR[n]:
										if y != k:
											if pos5[0] <= y[1] + 2*INSERT and  pos5[0] + vasouille >= y[1] and pos6[0] <= y[4] + 2*INSERT and  pos6[0] + vasouille >= y[4]:
												if pos1[1] <= y[0] + vasouille and pos6[1] <= pos4[0] + vasouille:
													outfile.write('reciprocal_translocation\tregion1_inv:\t'+n+'\t'+str(pos1[0])+'\t'+str(y[1])+'\tregion2_inv:\t'+n+'\t'+str(pos6[0])+'\t'+str(pos4[1])+'\n')
	#Now it's time to detect inversion of the first or the second fragments of reciprocal translocation
	if TYPE == '0':
		for n in dic_FF:
			for j in dic_FF[n]:
				pos1 = [j[0], j[1]]
				pos3 = [j[3], j[4]]
				pos2 = []
				pos4 = []
				if n in dic_FR:
					for k in dic_FR[n]:
						if pos3[0] <= k[4] + 2*INSERT and  pos3[0] + vasouille >= k[4]:
							pos2 = [k[0], k[1]]
							pos4 = [k[3], k[4]]
							if n in dic_RR:
								for x in dic_RR[n]:
									if pos2[0] + vasouille >= x[1] and pos2[0] <= x[1] + 2*INSERT:
										pos5 = [x[0], x[1]]
										pos6 = [x[3], x[4]]
										if n in dic_DEL:
											for y in dic_DEL[n]:
												if y[3] + vasouille >= pos6[1] and y[3] <= pos6[1] + 2*INSERT and y[1] <= pos1[0] + vasouille and y[1] + 2*INSERT >= pos1[0]:
													pos7 = [y[0], y[1]]
													pos8 = [y[3], y[4]]
													if pos2[0] + vasouille >= pos1[1]:#first fragment inversed
														outfile.write('reciprocal_translocation\tregion1_inv:\t'+n+'\t'+str(pos1[0])+'\t'+str(pos5[1])+'\tregion2:\t'+n+'\t'+str(pos8[0])+'\t'+str(pos4[1])+'\n')
													elif pos2[1] <= pos1[0] + vasouille:#second fragment inversed
														outfile.write('reciprocal_translocation\tregion1:\t'+n+'\t'+str(pos2[0])+'\t'+str(pos7[1])+'\tregion2_inv:\t'+n+'\t'+str(pos3[0])+'\t'+str(pos6[1])+'\n')
	##################################################################
	#####This part only works for inter chromosomal rearrangments#####
	##################################################################
	#Now it's time to detect simple translocations
	if TYPE == '0':
		for n in dic_DEL:
			for j in dic_DEL[n]:
				pos1 = [j[0], j[1]]
				pos2 = [j[3], j[4]]
				pos3 = []
				pos4 = []
				pos5 = []
				pos6 = []
				if n in dic_CHR_fr_rev:
					for k in dic_CHR_fr_rev[n]:
						if k[1] <= pos2[0] + vasouille and k[1] + 2*INSERT >= pos2[0]:
							pos4 = [k[0], k[1]]
							pos3 = [k[3], k[4]]
							if n in dic_CHR_rf_rev:
								for m in dic_CHR_rf_rev[n]:
									if m[0] + vasouille >= pos1[1] and m[0] <= pos1[1] + 2*INSERT and m[2] == k[2]:
										pos6 = [m[0],m[1]]
										pos5 = [m[3],m[4]]
										if pos5[1] <= pos3[0] + vasouille and pos5[1] + 2*INSERT >= pos3[0]:
											outfile.write('translocation\tregion:\t'+n+'\t'+str(pos6[0])+'\t'+str(pos4[1])+'\ttarget:\t'+k[2]+'\t'+str(pos5[1])+'\t'+str(pos3[0])+'\n')
				if n in dic_CHR_fr:
					for k in dic_CHR_fr[n]:
						if k[0] + vasouille >= pos1[1] and k[0] <= pos1[1] + 2*INSERT:
							pos4 = [k[0], k[1]]#inversion pos4 et pos3
							pos3 = [k[3], k[4]]
							if n in dic_CHR_rf:
								for m in dic_CHR_rf[n]:
									if m[1] <= pos2[0] + vasouille and m[1] + 2*INSERT >= pos2[0] and m[2] == k[2]:
										pos6 = [m[0],m[1]]#inversion pos6 et pos5
										pos5 = [m[3],m[4]]
										if pos3[1] <= pos5[0] + vasouille and pos3[1] + 2*INSERT >= pos5[0]:
											outfile.write('translocation\tregion:\t'+n+'\t'+str(pos4[0])+'\t'+str(pos6[1])+'\ttarget:\t'+k[2]+'\t'+str(pos3[1])+'\t'+str(pos5[0])+'\n')
	#Now it's time to detect simple duplications
	if TYPE == '6':
		for n in dic_CHR_rf:
			for j in dic_CHR_rf[n]:
				pos1 = [j[0], j[1]]
				pos2 = [j[3], j[4]]
				pos3 = []
				pos4 = []
				if n in dic_CHR_fr:
					for k in dic_CHR_fr[n]:
						if k[0] + vasouille >= pos1[1] and k[0] <= pos1[1] + 2*INSERT and j[2] == k[2]:
							pos3 = [k[0], k[1]]
							pos4 = [k[3], k[4]]
							if couv_med_reg(pos2[0], pos4[1], j[2]) >= (EXP_COV + float(EXP_COV)*PLOID):
								outfile.write('duplication\tregion:\t'+k[2]+'\t'+str(pos2[0])+'\t'+str(pos4[1])+'\ttarget:\t'+n+'\t'+str(pos1[1])+'\t'+str(pos3[0])+'\n')
	if TYPE == '7':
		for n in dic_CHR_fr_rev:
			for j in dic_CHR_fr_rev[n]:
				pos1 = [j[0], j[1]]
				pos2 = [j[3], j[4]]
				pos3 = []
				pos4 = []
				if n in dic_CHR_rf_rev:
					for k in dic_CHR_rf_rev[n]:
						if k[0] + vasouille >= pos1[1] and k[0] <= pos1[1] + 2*INSERT and j[2] == k[2]:
							pos3 = [k[0], k[1]]
							pos4 = [k[3], k[4]]
							if couv_med_reg(pos2[0], pos4[1], j[2]) >= (EXP_COV + float(EXP_COV)*PLOID):
								outfile.write('duplication\tregion:\t'+k[2]+'\t'+str(pos2[0])+'\t'+str(pos4[1])+'\ttarget:\t'+n+'\t'+str(pos1[1])+'\t'+str(pos3[0])+'\n')
	#Now it's time to detect simple translocations with inversion
	if TYPE == '0':
		for n in dic_DEL:
			for j in dic_DEL[n]:
				pos1 = [j[0], j[1]]
				pos2 = [j[3], j[4]]
				pos3 = []
				pos4 = []
				pos5 = []
				pos6 = []
				if n in dic_CHR_rr_rev:
					for k in dic_CHR_rr_rev[n]:
						if k[1] <= pos2[0] + vasouille and k[1] + 2*INSERT >= pos2[0]:
							pos4 = [k[0], k[1]]
							pos3 = [k[3], k[4]]
							if n in dic_CHR_ff_rev:
								for m in dic_CHR_ff_rev[n]:
									if m[0] + vasouille >= pos1[1] and m[0] <= pos1[1] + 2*INSERT and m[2] == k[2]:
										pos6 = [m[0],m[1]]
										pos5 = [m[3],m[4]]
										if pos3[1] <= pos5[0] + vasouille and pos3[1] + 2*INSERT >= pos5[0]:
											outfile.write('translocation\tregion_inv:\t'+n+'\t'+str(pos6[0])+'\t'+str(pos4[1])+'\ttarget:\t'+k[2]+'\t'+str(pos3[1])+'\t'+str(pos5[0])+'\n')
				if n in dic_CHR_ff:
					for k in dic_CHR_ff[n]:
						if k[0] + vasouille >= pos1[1] and k[0] <= pos1[1] + 2*INSERT:
							pos3 = [k[0], k[1]]
							pos4 = [k[3], k[4]]
							if n in dic_CHR_rr:
								for m in dic_CHR_rr[n]:
									if m[1] <= pos2[0] + vasouille and m[1] + 2*INSERT >= pos2[0] and m[2] == k[2]:
										pos5 = [m[0],m[1]]
										pos6 = [m[3],m[4]]
										if pos6[1] <= pos4[0] + vasouille and pos6[1] + 2*INSERT >= pos4[0]:
											outfile.write('translocation\tregion_inv:\t'+n+'\t'+str(pos3[0])+'\t'+str(pos5[1])+'\ttarget:\t'+k[2]+'\t'+str(pos6[1])+'\t'+str(pos4[0])+'\n')
	#Now it's time to detect simple duplications with inversion
	if TYPE == '8':
		for n in dic_CHR_rr:
			for j in dic_CHR_rr[n]:
				pos1 = [j[0], j[1]]
				pos2 = [j[3], j[4]]
				pos3 = []
				pos4 = []
				if n in dic_CHR_ff:
					for k in dic_CHR_ff[n]:
						if k[0] + vasouille >= pos1[1] and k[0] <= pos1[1] + 2*INSERT and j[2] == k[2]:
							pos3 = [k[0], k[1]]
							pos4 = [k[3], k[4]]
							if couv_med_reg(pos4[0], pos2[1], j[2]) >= (EXP_COV + float(EXP_COV)*PLOID):
								outfile.write('duplication\tregion_inv:\t'+k[2]+'\t'+str(pos4[0])+'\t'+str(pos2[1])+'\ttarget:\t'+n+'\t'+str(pos1[1])+'\t'+str(pos3[0])+'\n')
	if TYPE == '9':
		for n in dic_CHR_rr_rev:
			for j in dic_CHR_rr_rev[n]:
				pos1 = [j[0], j[1]]
				pos2 = [j[3], j[4]]
				pos3 = []
				pos4 = []
				if n in dic_CHR_ff_rev:
					for k in dic_CHR_ff_rev[n]:
						if k[0] + vasouille >= pos1[1] and k[0] <= pos1[1] + 2*INSERT and j[2] == k[2]:
							pos3 = [k[0], k[1]]
							pos4 = [k[3], k[4]]
							if couv_med_reg(pos4[0], pos2[1], j[2]) >= (EXP_COV + float(EXP_COV)*PLOID):
								outfile.write('duplication\tregion_inv:\t'+k[2]+'\t'+str(pos4[0])+'\t'+str(pos2[1])+'\ttarget:\t'+n+'\t'+str(pos1[1])+'\t'+str(pos3[0])+'\n')
	#Now it's time to detect simple reciprocal translocations
	if TYPE == '0':
		for n in dic_CHR_fr:
			dic_CHR_fr_rf = []
			for j in dic_CHR_fr[n]:
				pos1 = [j[0], j[1]]
				pos2 = [j[3], j[4]]
				if n in dic_CHR_rf:
					for k in dic_CHR_rf[n]:
						if pos1[0] <= k[1] + 2*INSERT and pos1[0] + vasouille >= k[1] and k[3] <= pos2[1] + 2*INSERT and k[3] + vasouille >= pos2[1] and j[2] == k[2]:
							pos3 = [k[0], k[1]]
							pos4 = [k[3], k[4]]
							dic_CHR_fr_rf.append([pos3[0], pos1[1], pos2[0], pos4[1], pos1[0], pos4[0], pos3[1], pos2[1], j[2]])
			L_done = []
			for x in dic_CHR_fr_rf:
				L_done.append(x)
				for y in dic_CHR_fr_rf:
					if not(y in L_done):
						if x[1] <= y[0] + vasouille and  x[3] <= y[2] + vasouille and x[8] == y[8]:
							outfile.write('reciprocal_translocation\tregion1:\t'+n+'\t'+str(x[4])+'\t'+str(y[6])+'\tregion2:\t'+x[8]+'\t'+str(x[5])+'\t'+str(y[7])+'\n')
						elif y[1] <= x[0] + vasouille and y[3] <= x[2] + vasouille and x[8] == y[8]:
							outfile.write('reciprocal_translocation\tregion1:\t'+n+'\t'+str(y[4])+'\t'+str(x[6])+'\tregion2:\t'+x[8]+'\t'+str(y[5])+'\t'+str(x[7])+'\n')
	#Now it's time to detect inversion of both fragments of reciprocal translocation
	if TYPE == '0':
		for n in dic_CHR_ff:
			for j in dic_CHR_ff[n]:
				pos1 = [j[0], j[1]]
				pos3 = [j[3], j[4]]
				pos2 = []
				pos4 = []
				if n in dic_CHR_rr:
					for k in dic_CHR_rr[n]:
						if pos1[0] <= k[1] + 2*INSERT and  pos1[0] + vasouille >= k[1] and pos3[0] <= k[4] + 2*INSERT and  pos3[0] + vasouille >= k[4] and k[2] == j[2]:
							pos2 = [k[0], k[1]]
							pos4 = [k[3], k[4]]
							for x in dic_CHR_ff[n]:
								if x != j:
									pos5 = [x[0], x[1]]
									pos6 = [x[3], x[4]]
									for y in dic_CHR_rr[n]:
										if y != k:
											if pos5[0] <= y[1] + 2*INSERT and  pos5[0] + vasouille >= y[1] and pos6[0] <= y[4] + 2*INSERT and  pos6[0] + vasouille >= y[4] and x[2] == y[2] and x[2] == k[2]:
												if pos1[1] <= y[0] + vasouille and pos6[1] <= pos4[0] + vasouille:
													outfile.write('reciprocal_translocation\tregion1_inv:\t'+n+'\t'+str(pos1[0])+'\t'+str(y[1])+'\tregion2_inv:\t'+k[2]+'\t'+str(pos6[0])+'\t'+str(pos4[1])+'\n')
	#Now it's time to detect inversion of the first or the second fragments of reciprocal translocation
	if TYPE == '0':
		for n in dic_CHR_ff:
			for j in dic_CHR_ff[n]:
				pos1 = [j[0], j[1]]
				pos3 = [j[3], j[4]]
				pos2 = []
				pos4 = []
				if j[2] in dic_CHR_fr_rev:
					for k in dic_CHR_fr_rev[j[2]]:
						if pos3[0] <= k[1] + 2*INSERT and  pos3[0] + vasouille >= k[1] and n == k[2]:
							pos4 = [k[0], k[1]]
							pos2 = [k[3], k[4]]
							if n in dic_CHR_rr:
								for x in dic_CHR_rr[n]:
									if pos2[0] + vasouille >= x[1] and pos2[0] <= x[1] + 2*INSERT and j[2] == x[2]:
										pos5 = [x[0], x[1]]
										pos6 = [x[3], x[4]]
										if j[2] in dic_CHR_rf_rev:
											for y in dic_CHR_rf_rev[j[2]]:
												if y[0] + vasouille >= pos6[1] and y[0] <= pos6[1] + 2*INSERT and y[4] <= pos1[0] + vasouille and y[4] + 2*INSERT >= pos1[0] and n == y[2]:
													pos8 = [y[0], y[1]]
													pos7 = [y[3], y[4]]
													if pos2[0] + vasouille >= pos1[1]:#first fragment inversed
														outfile.write('reciprocal_translocation\tregion1_inv:\t'+n+'\t'+str(pos1[0])+'\t'+str(pos5[1])+'\tregion2:\t'+j[2]+'\t'+str(pos8[0])+'\t'+str(pos4[1])+'\n')
													elif pos2[1] <= pos1[0] + vasouille:#second fragment inversed
														outfile.write('reciprocal_translocation\tregion1:\t'+n+'\t'+str(pos2[0])+'\t'+str(pos7[1])+'\tregion2_inv:\t'+j[2]+'\t'+str(pos3[0])+'\t'+str(pos6[1])+'\n')
	outfile.close()


def couv_med_reg(deb, fin, chr):
	"""
		:param chr: The chromosome name
		:type chr: str
		:param deb: the start position
		:type deb: int
		:param end: The end position
		:type end: int
	"""
	subListNoGap = []
	for site in covChr[chr][deb-1:fin]:
		if site:
			subListNoGap.append(site)
	if len(subListNoGap) == 0:
		return 0
	else:
		return mediane(subListNoGap)
		

def couv_med_reg_del(deb, fin, chr):
	"""
		:param chr: The chromosome name
		:type chr: str
		:param deb: the start position
		:type deb: int
		:param end: The end position
		:type end: int
	"""
	subListNoGap = []
	for site in covChr[chr][deb-1:fin]:
		if site:
			subListNoGap.append(site)
	if len(subListNoGap) <= 0.5*(fin-deb):
		return 0
	else:
		return mediane(subListNoGap)
		

def mediane(LISTE):
	L = sorted(LISTE)
	N = len(L)
	n = N/2.0
	p = int(n)
	if n == 1:
		return (L[0])
	elif n == p:
		return (L[p-1]+L[p])/2.0
	else:
		return L[p]


def main(job):
	indent_discord(job[0]['ff'], job[0]['frf'], job[0]['rr'], job[0]['ins'], job[0]['delet'], job[0]['chr_rr'], job[0]['chr_fr'], job[0]['chr_rf'], job[0]['chr_ff'], job[0]['insert'], job[1], job[0]['medCoverage'], job[0]['ploid'], job[2])
	return 0

def worker(job):
	print job
	try:
		codeError = main(job)
	except Exception as e:
		print e
		codeError = 1
	finally:
		return codeError
		
		
def __main__():
	#Parse Command Line
	parser = optparse.OptionParser(usage="python %prog [options]\n\nProgram designed by Guillaume MARTIN : guillaume.martin@cirad.fr\n\n"
	"This script try to identify signature for structural variation (this is a wrapper for idend_sv.py)")
	# Wrapper options.
	parser.add_option( '', '--frf', dest='frf', default='not_filled', help='The discordant fr or rf file depending on expected orientation')
	parser.add_option( '', '--ff', dest='ff', default='not_filled', help='The discordant ff file')
	parser.add_option( '', '--rr', dest='rr', default='not_filled', help='The discordant rr file')
	parser.add_option( '', '--ins', dest='ins', default='not_filled', help='The discordant ins file')
	parser.add_option( '', '--delet', dest='delet', default='not_filled', help='The discordant del file')
	parser.add_option( '', '--chr_rr', dest='chr_rr', default='not_filled', help='The discordant chr_rr file')
	parser.add_option( '', '--chr_fr', dest='chr_fr', default='not_filled', help='The discordant chr_fr file')
	parser.add_option( '', '--chr_rf', dest='chr_rf', default='not_filled', help='The discordant chr_rf file')
	parser.add_option( '', '--chr_ff', dest='chr_ff', default='not_filled', help='The discordant chr_ff file')
	parser.add_option( '', '--chr', dest='chr', default='not_filled', help='The tabulated file containing in col 1 : chromosome name, col 2: chromosome length. A line for each chromosomes')
	parser.add_option( '', '--covf', dest='covf', default='not_filled', help='The coverage file calculated in "5_calc_stat"')
	parser.add_option( '', '--orient', dest='orient', default='rf', help='The expected orientation: rf or fr, [default: %default]')
	parser.add_option( '', '--insert', dest='insert', default=5000, help='The expected insert size, [default: %default]')
	parser.add_option( '', '--exp_cov', dest='exp_cov', default='not_filled', help='The expected coverage (float)')
	parser.add_option( '', '--ploid', dest='ploid', default=0.33, help='Multiplicator for coverage variation detection in SV identification (ex : If homozygous duplication expected in diploid: expected = coverage + coverage*1, if heterozygous duplication expected in diploid: expected = coverage + coverage*0.5). Choose a value lower than the expected one')
	parser.add_option( '', '--type', dest='type', default='0123456789', help='The type of SV searched: q')
	parser.add_option( '', '--thread', dest='thread', default='1', help='The thread number used for search (integer), [default: %default]')
	parser.add_option( '', '--out', dest='out', default='SV_detected.tab', help='Output file')
	parser.add_option( '', '--config', dest='config', default=None)
	(options, args) = parser.parse_args()
	
	proc = int(options.thread)
	
	pathname = os.path.dirname(sys.argv[0])
	loca_programs = ConfigParser.RawConfigParser()
	loca_programs.read(pathname+'/loca_programs.conf')
	
	if options.frf == 'not_filled':
		sys.exit('Please provide an argument for --frf')
	if options.ff == 'not_filled':
		sys.exit('Please provide an argument for --ff')
	if options.rr == 'not_filled':
		sys.exit('Please provide an argument for --rr')
	if options.chr_rr == 'not_filled':
		sys.exit('Please provide an argument for --chr_rr')
	if options.chr_rf == 'not_filled':
		sys.exit('Please provide an argument for --chr_rf')
	if options.chr_fr == 'not_filled':
		sys.exit('Please provide an argument for --chr_fr')
	if options.chr_ff == 'not_filled':
		sys.exit('Please provide an argument for --chr_ff')
	if options.ins == 'not_filled':
		sys.exit('Please provide an argument for --ins')
	if options.chr_rr == 'not_filled':
		sys.exit('Please provide an argument for --chr_rr')
	if options.out == 'not_filled':
		sys.exit('Please provide an argument for --out')
		
	if not options.config:
		if options.exp_cov == 'not_filled':
			sys.exit('Please provide an argument for --exp_cov')
		if options.covf == 'not_filled':
			sys.exit('Please provide an argument for --covf')
		
	args = {}
	if options.orient == 'rf':
		args['ff'] = options.ff
		args['rr'] = options.rr
	elif options.orient == 'fr':
		args['ff'] = options.rr 
		args['rr'] = options.ff
	else:
		raise ValueError("Unrecognised orientation")
		
	if options.config:
		config = ConfigParser.RawConfigParser()
		config.read(options.config)
		args['insert'] = config.getfloat('Calc_coverage', 'median_insert')
		args['medCoverage'] = config.getfloat('Calc_coverage', 'median_coverage')
		args['covFile'] = config.get('Calc_coverage', 'out')
		args['ploid'] = config.getfloat('General','ploid')
		args['chrFile'] = config.get('General', 'chr')
	else:
		args['insert'] = float(options.insert)
		args['medCoverage'] = float(options.exp_cov)
		args['covFile'] = options.covf
		args['ploid'] = float(options.ploid)
		args['chrFile'] = options.chr
	args['frf'] = options.frf
	args['chr_rr'] = options.chr_rr
	args['chr_rf'] = options.chr_rf
	args['chr_fr'] = options.chr_fr
	args['chr_ff'] = options.chr_ff
	args['ins'] = options.ins
	args['delet'] = options.delet
	
	
	listJobs = []
	oufileNames = []
	for i in range(len(options.type)):
		outfile = "tmp_outfile_"+str(i)
		listJobs.append([args, outfile, options.type[i]])
		oufileNames.append(outfile)

	readCov(args['chrFile'], args['covFile'])
	
	pool = mp.Pool(processes=proc)
	resultsJobs = pool.map(worker, listJobs)
	
	for rslt in resultsJobs:
		print rslt
			

	with open(options.out, 'w') as outfile:		
		for fname in oufileNames:
			if os.path.isfile(fname):
				with open(fname) as infile:
					for line in infile:
						outfile.write(line)
				os.remove(fname)
			else:
				print("File "+fname+" not found.")

	if options.config:
		config = ConfigParser.RawConfigParser()
		config.read(options.config)
		config.set('Ident_discord','frf', options.frf)
		config.set('Ident_discord','ff', options.ff)
		config.set('Ident_discord','rr', options.rr)
		config.set('Ident_discord','ins', options.ins)
		config.set('Ident_discord','delet', options.delet)
		config.set('Ident_discord','chr_rr', options.chr_rr)
		config.set('Ident_discord','chr_fr', options.chr_fr)
		config.set('Ident_discord','chr_rf', options.chr_rf)
		config.set('Ident_discord','chr_ff', options.chr_ff)
		with open(options.config, 'wb') as configfile:
			config.write(configfile)		
	
	sys.exit()
	i = 0
	liste_tmp = []
	liste_id = []
	listeJobs = []
	# liste_job = []
	while i < 10:
		temp = options.out+'_'+str(i)
		liste_tmp.append(temp)
		# print temp
		if options.config:
			# liste_id.append("%s %s/ident_SV.py --frf %s --ff %s --rr %s --ins %s --delet %s --chr_rr %s --chr_fr %s --chr_rf %s --chr_ff %s --chr %s --covf %s --orient %s --insert %s --exp_cov %s --ploid %s --out %s --config %s --type %s" % (loca_programs.get('Programs','python'), ScriptPath, options.frf, options.ff, options.rr, options.ins, options.delet, options.chr_rr, options.chr_fr, options.chr_rf, options.chr_ff, options.chr, options.covf, options.orient, options.insert, options.exp_cov, options.ploid, temp, options.config, str(i)))
			args = [ScriptPath, options.frf, options.ff, options.rr, options.ins, options.delet, options.chr_rr, options.chr_fr, options.chr_rf, options.chr_ff, options.chr, options.covf, options.orient, options.insert, options.exp_cov, options.ploid, temp, options.config, str(i)]
			type = str(i)
			listeJobs.append(args)
		else:
			# liste_id.append("%s %s/ident_SV.py --frf %s --ff %s --rr %s --ins %s --delet %s --chr_rr %s --chr_fr %s --chr_rf %s --chr_ff %s --chr %s --covf %s --orient %s --insert %s --exp_cov %s --ploid %s --out %s --type %s" % (loca_programs.get('Programs','python'), ScriptPath, options.frf, options.ff, options.rr, options.ins, options.delet, options.chr_rr, options.chr_fr, options.chr_rf, options.chr_ff, options.chr, options.covf, options.orient, options.insert, options.exp_cov, options.ploid, temp, str(i)))
			args = [ScriptPath, options.frf, options.ff, options.rr, options.ins, options.delet, options.chr_rr, options.chr_fr, options.chr_rf, options.chr_ff, options.chr, options.covf, options.orient, options.insert, options.exp_cov, options.ploid, temp, str(i)]
			listeJobs.append(args)
		i += 1
	
	liste_process = []
	
	pool = multiprocessing.Pool(processes=proc)
	resultsJobs = pool.map(worker, liste_id)
	
	# for n in liste_id:
		# t = multiprocessing.Process(target=run_job, args=(n, 'Bug lauching indent_SV.py',))
		# liste_process.append(t)
		# if len(liste_process) == proc:
			# # Starts threads
			# for process in liste_process:
				# process.start()
			# # This blocks the calling thread until the thread whose join() method is called is terminated.
			# for process in liste_process:
				# process.join()
			# #the processes are done
			# liste_process = []
	# if liste_process:
		# # Starts threads
		# for process in liste_process:
			# process.start()
		# # This blocks the calling thread until the thread whose join() method is called is terminated.
		# for process in liste_process:
			# process.join()
		# #the processes are done
		# liste_process = []
	
	if options.config:
		config = ConfigParser.RawConfigParser()
		config.read(options.config)
		config.set('Ident_discord','frf', options.frf)
		config.set('Ident_discord','ff', options.ff)
		config.set('Ident_discord','rr', options.rr)
		config.set('Ident_discord','ins', options.ins)
		config.set('Ident_discord','delet', options.delet)
		config.set('Ident_discord','chr_rr', options.chr_rr)
		config.set('Ident_discord','chr_fr', options.chr_fr)
		config.set('Ident_discord','chr_rf', options.chr_rf)
		config.set('Ident_discord','chr_ff', options.chr_ff)
		with open(options.config, 'wb') as configfile:
			config.write(configfile)
	# for n in liste_job:
		# cherche_error('IDENT_SV.o'+n)
		# os.system('rm IDENT_SV.o'+n)
	mot = 'cat '
	for n in liste_tmp:
		mot = mot + n + ' '
	mot = mot + '> ' + options.out
	os.system(mot)
	for n in liste_tmp:
		os.remove(n)
		
if __name__ == "__main__": __main__()