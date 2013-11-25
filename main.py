#!/usr/bin/env python
# -*- coding: latin-1 -*-
'''Function to evaluate the activity of each species during the simulation, 
   catalyst substrate product or nothing
'''

import sys, os # Standard library
import datetime as dt
import linecache as lc
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

#ÊInput parameters definition 
if __name__ == '__main__':
	parser = ArgumentParser(
				description='Main script of ACS analysis.'
				, epilog='''ACS ANALYSIS Main File. ''') 
	parser.add_argument('-i', '--initanal', type=int, help='Analysis of the initial structures', choices=[0,1], default=0)
	parser.add_argument('-e', '--exhaustive', type=int, help='Analysis of the simulation', choices=[0,1], default=1)
	parser.add_argument('-d', '--decay', type=int, help='Decay time', default=0)
	parser.add_argument('-m', '--maxDim', help='Max Dimension of the system', default='4', type=int)
	parser.add_argument('-p', '--strPath', help='Path where files are stored', default='./')
	parser.add_argument('-r', '--resFolder', help='Name of the result folder', default='res')
	
	args = parser.parse_args()
	
	# Create absolute paths
	strPath = os.path.abspath(args.strPath)
	tmpDirs = sort(os.listdir(strPath))
	
	# Goes into the simulation folder
	os.chdir(strPath)
	
	# Create stas folders
	ndn = '_0_new_allStatResults'
	newdirAllResults = os.path.join(strPath, ndn)
	if not os.path.isdir(newdirAllResults):
		try: os.mkdir(newdirAllResults)
		except: print "Impossible to create statistic directory", newdirAllResults; sys.exit(1)
	print "|- Simulation Folder: ", strPath
	# System Type
	_CLOSE_ = 0
	_PROTO_ = 1
	_CSTR_ = 2
	
	if args.initanal == 1:
		fname_initRafRes = os.path.join(newdirAllResults, '0_initRafAnalysis.csv')
		fname_initRafResLIST = os.path.join(newdirAllResults, '0_initRafAnalysisLIST.csv')
		fname_initRafResALL = os.path.join(newdirAllResults, '0_initRafAnalysisALL.csv')
		fid_initRafRes = open(fname_initRafRes, 'w')
		fid_initRafResLIST = open(fname_initRafResLIST, 'w')
		fid_initRafResALL = open(fname_initRafResALL, 'w')
		strToWrite = "Folder\tP\tM\tRAFsize\tClosure\tCats\tuRAF\n"
		fid_initRafRes.write(strToWrite)
	
	for tmpDir in tmpDirs:
		
		totDirName = os.path.join(strPath,tmpDir)
		if os.path.isdir(totDirName):
			# Move to the directory 
			os.chdir(totDirName)
			
			# SIMULATION DYNAMICS FOLDER
			resDirPath = os.path.abspath(os.path.join("./", args.resFolder))
			print " \- Results Folder: ", totDirName
			if os.path.isdir(resDirPath):
				
				# Analysis of the initial structures 
				conf = readfiles.readConfFile(totDirName) # Configuration file upload
				# System type recognition
				if (conf[6] == 0) & (conf[7] > 0): sysType = _PROTO_
				elif (conf[6] > 0) & (conf[7] == 0): sysType = _CSTR_
				elif (conf[6] == 0) & (conf[7] == 0): sysType = _CLOSE_
				foodList = dm.generateFluxList(totDirName, sysType)
								
				# Initial Analysis are turned ON (RAF ANALYSIS)
				if args.initanal == 1:
					print "   |- LOADING init structures..."
					rcts = readfiles.loadAllData(totDirName,'_acsreactions.csv') # reaction file upload
					cats = readfiles.loadAllData(totDirName,'_acscatalysis.csv') # catalysis file upload
					# REAL RAF COMPUTATION (!!!!!! 1 is average connectivity)
					raf.rafComputation(fid_initRafRes, fid_initRafResALL, fid_initRafResLIST, tmpDir, conf[9], 1, rcts, cats, foodList, args.maxDim)
					
				os.chdir(resDirPath)
				print "  \-Folder ", resDirPath
				# Find the number of generations
				numberOfGen = len(glob.glob(os.path.join(resDirPath,'times_*')))
				# For each generation
				os.chdir(resDirPath)
				
				# Analysis of the dynamics
				if args.exhaustive == 1:
					for ngen in range(1,numberOfGen+1):
						print "\t\t\tGeneration ", ngen
						strZeros = readfiles.zeroBeforeStrNum(ngen, numberOfGen)
						# Saving File Creation
						if args.exhaustive == 1:
							fName = 'RAF_dynamics_analysis_gen_' + strZeros + str(ngen) + '.csv'
							fname_dynRafRes = os.path.join(newdirAllResults, fName)
							fid_initRafRes = open(fname_dynRafRes, 'w')
							strToWrite = "t\t#CL\t#RAF\tRAF\n"
							fid_initRafRes.write(strToWrite)
						# reactions file loading
						if ngen == 1:
					  		strRctZero = 'reactions_' + strZeros + str(0) + '*'
					  		strCatZero = 'catalysis_' + strZeros + str(0) + '*'
					  		rctFilesZero = sorted(glob.glob(os.path.join(resDirPath,strRctZero)))
					  		catFilesZero = sorted(glob.glob(os.path.join(resDirPath,strCatZero)))
						strRct = 'reactions_' + strZeros + str(ngen) + '*'  
						strCat = 'catalysis_' + strZeros + str(ngen) + '*'  
					  	# Searching for files
					  	rctFiles = sorted(glob.glob(os.path.join(resDirPath,strRct)))
					  	catFiles = sorted(glob.glob(os.path.join(resDirPath,strCat)))
					  	
					  	if ngen == 1: # If the generation is the first one the 0000 files must be uploaded too
				  			rctFiles = rctFilesZero + rctFiles
				  			catFiles = catFilesZero + catFiles
				  		
				  		sngTime = conf[2] / (float(len(rctFiles)) - 1)
					  	actTime = 0
				  		
				  		for id, sngFile in enumerate(rctFiles):
				  			# Load files
				  			rcts = readfiles.loadAllData(totDirName,sngFile) # reaction file upload
							cats = readfiles.loadAllData(totDirName,catFiles[id]) # catalysis file upload
							foodList = dm.generateFluxList(totDirName, sysType)
							# Reshape rcts according to the reactions really occurred
							procrcts = rcts[rcts[:,5] > 0,:]
							proccats = cats[cats[:,3] > 0,:]
							R = raf.rafDynamicComputation(fid_initRafRes, actTime, procrcts[:,0:5], proccats[:,0:5], foodList)
							#print R
							actTime += sngTime
						
						# If decay analysis is TURNED ON (rections are added step by step to the graph and they are removed 
						# if they do not occur twice within the decay time
						
					
					
	if args.initanal == 1: 
		fid_initRafRes.close()
		fid_initRafResLIST.close()
		fid_initRafResALL.close()
					
				  	
				  	


