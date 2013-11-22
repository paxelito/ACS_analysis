#!/usr/bin/env python
# -*- coding: latin-1 -*-
'''File to convert Wim files in my files. 
'''

import sys, os # Standard library
import datetime as dt
import linecache as lc
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
	parser.add_argument('-p', '--strPath', help='Path where files are stored', default='./')
	parser.add_argument('-w', '--wimFile', help='Wim File', default='RAF.crs')
	parser.add_argument('-o', '--strOut', help='Path for the file storing', default='./')
	
	args = parser.parse_args()

	# Create stas folders
	ndn = '_0_new_allStatResults'
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
			except:
				procedure = 0
				print foodList
		if (procedure == 3) & (line[0] != '<'):
			l = line.split()
			if len(l) > 0:
				# CONDENSATION
				if l[2] == '+': 
					if rct == 0: rcts = np.array([[int(rct), int(0), int(speciesList.index(l[5])), int(speciesList.index(l[1])), int(speciesList.index(l[3])), int(0), int(0), int(0)]])
					else: rcts = np.vstack([rcts,(int(rct), int(0), int(speciesList.index(l[5])), int(speciesList.index(l[1])), int(speciesList.index(l[3])), int(0), int(0), int(0))])
				if l[2] == '->': # CLEAVAGE
					if rct == 0: rcts = np.array([[int(rct), int(1), int(speciesList.index(l[1])), int(speciesList.index(l[3])), int(speciesList.index(l[5])), int(0), int(0), int(0)]])
					else: rcts = np.vstack([rcts,(int(rct), int(1), int(speciesList.index(l[1])), int(speciesList.index(l[3])), int(speciesList.index(l[5])), int(0), int(0), int(0))])
				
				# CATALYST
				catNums = len(l) - 7
				if catNums == 1:
					if rct == 0: cats = np.array([[int(rct), int(speciesList.index(l[6][1:len(l[6])-1])), int(rct), int(0), int(5), int(25), int(5), int(1)]])
					else: cats = np.vstack([cats,(int(rct), int(speciesList.index(l[6][1:len(l[6])-1])), int(rct), int(0), int(5), int(25), int(5), int(1))])					
				if catNums > 1:
					for i in range(catNums):
						if i == 0:
							if rct == 0: cats = np.array([[int(rct), int(speciesList.index(l[6][1:len(l[6])])), int(rct), int(0), int(5), int(25), int(5), int(1)]])
							else: cats = np.vstack([cats,(int(rct), int(speciesList.index(l[6][1:len(l[6])])), int(rct), int(0), int(5), int(25), int(5), int(1))])							
						elif i == (catNums-1):
							cats = np.vstack([cats,(int(rct), int(speciesList.index(l[6+catNums-1][0:len(l[6+catNums-1])-1])), int(rct), int(0), int(5), int(25), int(5), int(1))])
						else: 
							cats = np.vstack([cats,(int(rct), int(speciesList.index(l[6+i])), int(rct), int(0), int(5), int(25), int(5), int(1))])
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
