#!/usr/bin/env python
# -*- coding: latin-1 -*-
'''File to convert Wim files in my files. 
'''

import sys, os # Standard library
import datetime as dt
import linecache as lc
import random as ran
from copy import deepcopy
import glob
from argparse import ArgumentParser
import numpy as np # Scientific library
from numpy import * 

try:
    from pylab import *
except:
    pass
   
from lib.IO import * 
from lib.graph import raf
from lib.dyn import dynamics as dm


#æInput parameters definition 
if __name__ == '__main__':
	parser = ArgumentParser(
				description='From WIM format to my format'
				, epilog='''From WIM format to my format. ''') 
	parser.add_argument('-p', '--strPath', help='Path where files are stored (def: ./)', default='./')
	parser.add_argument('-w', '--wimFile', help='Wim File (def: RAF.crs)', default='RAF.crs') 
	parser.add_argument('-o', '--strOut', help='Path for the file storing (def: ./)', default='./')
	parser.add_argument('-k', '--rates', help='catalysis reaction rates (deafult: [50, 25, 50])	', nargs=3, metavar=('k_cond', 'k_clea', 'k_cpx'), default= [50, 25, 50])
	parser.add_argument('-s', '--sponstrate', help='spontaneous reactions rate (Def: 1e-6)', type=float, default=0)
	parser.add_argument('-t', '--sponstrateratio', help='spontaneous reactions ratio with respect to not spontanoues one (Def: 0', type=float, default=0)
	parser.add_argument('-f', '--framework', help='System Framework (0: close, 1: proto, 2: cstr - Def: 1', type=int, default= 1)
	parser.add_argument('-n', '--initialNumerosity', help='initial numerosity (Def: 600)', type=int, default=600)
	parser.add_argument('-v', '--volume', help='volume (Def: 1e-18)', type=float, default=1e-18)
	parser.add_argument('-r', '--resFolder', help='Folder of simulation outputs (def: res) ', default="res")
	parser.add_argument('-d', '--lastFiringDiskSpeciesID', help='Last ID species of the firing disk (def: 29) ', type=int, default=29)
	
	args = parser.parse_args()
	
	_AVOGADRO_ = 6.022141e23

	# Create stas folders
	ndn = '_0_Hordijk_raf'
	newdirAllResults = os.path.join(args.strOut, ndn)
	if not os.path.isdir(newdirAllResults):
		try:
			os.mkdir(newdirAllResults)
		except:
			print "Impossible to create statistic directory", newdirAllResults; sys.exit(1)
	
	print "|- Output Folder: ", args.strOut
	fname_initRafRes = os.path.join(newdirAllResults, '0_initRafAnalysis.csv')
	fname_initRafResLIST = os.path.join(newdirAllResults, '0_initRafAnalysisLIST.csv')
	fname_initRafResALL = os.path.join(newdirAllResults, '0_initRafAnalysisALL.csv')
	fid_initRafRes = open(fname_initRafRes, 'w')
	fid_initRafResLIST = open(fname_initRafResLIST, 'w')
	fid_initRafResALL = open(fname_initRafResALL, 'w')
	strToWrite = "Folder\tP\tAC\tM\tRAFsize\tClosure\tCats\tuRAF\n"
	fid_initRafRes.write(strToWrite)
	
	f = open(args.wimFile)
	lines = f.readlines()
	procedure = 0
	speciesList = []
	foodList = []
	rct = 0
	catID = 0
	
	for line in lines:
		if line[0] == '<':
			if line[0:11] == '<molecules>': 
				print '|- Molecules Processing'
				procedure = 1
			if line[0:10] == '<food set>': 
				print '|- Food Processing'
				procedure = 2
			if line[0:11] == '<reactions>': 
				print '|- Reactions Processing'
				procedure = 3
		# SPECIES
		if (procedure == 1) & (line[0] != '<'):
			l = line.split()
			try:
				speciesList.append(l[1])
			except:
				procedure = 0
				#print speciesList
		# FOODLIST
		if (procedure == 2) & (line[0] != '<'):
			l = line.split()
			try:
				foodList.append(speciesList.index(l[1]))
				food = deepcopy(foodList)
			except:
				procedure = 0
				print foodList
		if (procedure == 3) & (line[0] != '<'):
			l = line.split()
			if len(l) > 0:
				if ran.random() < args.sponstrateratio: kspont = args.sponstrate
				else: kspont = 0
				# CONDENSATION
				if l[2] == '+': # CONDENSATION
					if rct == 0: rcts = np.array([[int(rct), int(0), int(speciesList.index(l[5])), int(speciesList.index(l[1])), int(speciesList.index(l[3])), int(0), int(0), kspont]], dtype=float)
					else: rcts = np.vstack([rcts,(int(rct), int(0), int(speciesList.index(l[5])), int(speciesList.index(l[1])), int(speciesList.index(l[3])), int(0), int(0), kspont)])
				if l[2] == '->': # CLEAVAGE
					if rct == 0: rcts = np.array([[int(rct), int(1), int(speciesList.index(l[1])), int(speciesList.index(l[3])), int(speciesList.index(l[5])), int(0), int(0), kspont]], dtype=float)
					else: rcts = np.vstack([rcts,(int(rct), int(1), int(speciesList.index(l[1])), int(speciesList.index(l[3])), int(speciesList.index(l[5])), int(0), int(0), kspont)])
				
				# CATALYST
				catNums = len(l) - 7
				if catNums == 1:
					if rct == 0: cats = np.array([[int(catID), int(speciesList.index(l[6][1:len(l[6])-1])), int(rct), int(0), int(args.rates[0]), int(args.rates[1]), int(args.rates[0]), int(2)]])
					else: cats = np.vstack([cats,(int(catID), int(speciesList.index(l[6][1:len(l[6])-1])), int(rct), int(0), int(args.rates[0]), int(args.rates[1]), int(args.rates[0]), int(2))])					
					catID += 1
				if catNums > 1:
					for i in range(catNums):
						if i == 0:
							if rct == 0: cats = np.array([[int(catID), int(speciesList.index(l[6][1:len(l[6])])), int(rct), int(0), int(args.rates[0]), int(args.rates[1]), int(args.rates[0]), int(2)]])
							else: cats = np.vstack([cats,(int(catID), int(speciesList.index(l[6][1:len(l[6])])), int(rct), int(0), int(args.rates[0]), int(args.rates[1]), int(args.rates[0]), int(2))])							
						elif i == (catNums-1):
							cats = np.vstack([cats,(int(catID), int(speciesList.index(l[6+catNums-1][0:len(l[6+catNums-1])-1])), int(rct), int(0), int(args.rates[0]), int(args.rates[1]), int(args.rates[2]), int(1))])
						else: 
							cats = np.vstack([cats,(int(catID), int(speciesList.index(l[6+i])), int(rct), int(0), int(args.rates[0]), int(args.rates[1]), int(args.rates[2]), int(1))])
						catID += 1
				rct += 1
	
	print "|- RAF ANALYSIS... "
	rafsets = raf.rafComputation(fid_initRafRes, fid_initRafResALL, fid_initRafResLIST, 'tmpDir', 0, rct/float(len(speciesList)), rcts, cats, foodList, 10)
	print "|- CLOSURE DIMENSION: ", len(rafsets[0])
	print "|- RAF DIMENSION: ", len(rafsets[2])
	print "|- RAF UNIQUE DIMENSION: ", rafsets[4]
	f.close()
	
	# Close text files
	fid_initRafRes.close()
	fid_initRafResLIST.close()
	fid_initRafResALL.close()
	
	# CREATE DYNAMIC EXPERIMENTs FOLDER
	expsFolder = os.path.join(args.strOut, "exps")
	if not os.path.isdir(expsFolder):
		try: os.mkdir(expsFolder)
		except: print "Impossible to create exps directory", expsFolder; sys.exit(1)
		
	saveStrFolder = expsFolder
	
	fname_run = os.path.join(saveStrFolder, 'simulation.sh')
	fid_run = open(fname_run, 'w')	
	
	for cx in [29]:
		for sx in [0, 1e-3, 1]:
			if sx > 0:
				for rx in [1.0, 0.5]:
					# CREATE DYNAMIC EXPERIMENTs FOLDER
					folder = str(cx) + "_" + str(rx) + "_" + str(sx)
					s_expsFolder = os.path.join(saveStrFolder, folder)
					if not os.path.isdir(s_expsFolder):
						try: os.mkdir(s_expsFolder)
						except: print "Impossible to create exps directory", s_expsFolder; sys.exit(1)
						
					print "\t|- Save species file"
					fname_initSpecies = os.path.join(s_expsFolder, '_acsspecies.csv')
					fid_initSpecies = open(fname_initSpecies, 'w')	
					initConc = args.initialNumerosity / (args.volume * _AVOGADRO_)
					for id, species in enumerate(speciesList):
						lastc = "0"
						if id in food: lastc = "1"
						if id > cx: initConc = 0
						str2w = str(id) + "\t" + species+ "\t"  + str(initConc) + "\t1\t1\t0\t0\t1\t0\t0\t0\t0\t0\t0\t" + lastc + "\n"
						fid_initSpecies.write(str2w)   
					fid_initSpecies.close()
	
					writefiles.write_acsms_file(s_expsFolder, lfds=cx, sr=rx, K_spont_diss=sx/10000.,K_spont_ass=sx, alf='10', mrevrct=1)
					writefiles.writeAllFilesAndCreateResFolder(s_expsFolder, args.resFolder, cats, rcts, food, spontRatio=rx, kspontass=sx, kspontdiss=sx/10000.)
					str2w = "echo 'Simulation " + folder + "'\nnice ./carness ./" + folder + "/ " + folder + "/res/ " + folder + "/ > 0_log_" + folder + ".log\n"
					fid_run.write(str2w) 
			else:
				# CREATE DYNAMIC EXPERIMENTs FOLDER
				folder = str(cx) + "_" + str(0) + "_" + str(sx)
				s_expsFolder = os.path.join(saveStrFolder, folder)
				if not os.path.isdir(s_expsFolder):
					try: os.mkdir(s_expsFolder)
					except: print "Impossible to create exps directory", s_expsFolder; sys.exit(1)
					
				print "\t|- Save species file"
				fname_initSpecies = os.path.join(s_expsFolder, '_acsspecies.csv')
				fid_initSpecies = open(fname_initSpecies, 'w')	
				initConc = args.initialNumerosity / (args.volume * _AVOGADRO_)
				for id, species in enumerate(speciesList):
					lastc = "0"
					if id in food: lastc = "1"
					if id > cx: initConc = 0
					str2w = str(id) + "\t" + species+ "\t"  + str(initConc) + "\t1\t1\t0\t0\t1\t0\t0\t0\t0\t0\t0\t" + lastc + "\n"
					fid_initSpecies.write(str2w)   
				fid_initSpecies.close()

				writefiles.write_acsms_file(s_expsFolder, lfds=cx, sr=0, K_spont_diss=sx/10000.,K_spont_ass=sx, alf='10', mrevrct=1)
				writefiles.writeAllFilesAndCreateResFolder(s_expsFolder, args.resFolder, cats, rcts, food, spontRatio=0, kspontass=sx, kspontdiss=sx/10000.)
				str2w = "echo 'Simulation " + folder + "'\nnice ./carness ./" + folder + "/ " + folder + "/res/ " + folder + "/ > 0_log_" + folder + ".log\n"
				fid_run.write(str2w) 
					
	
	fid_run.close()
	
	# Save WIM structures
	# Create species file
	
