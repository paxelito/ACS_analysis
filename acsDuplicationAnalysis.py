#!/usr/bin/env python
# -*- coding: latin-1 -*-
'''Script to compute the successive cell division times and the value of each molecule in proximity of the cell division
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

def zeroBeforeStrNum(tmpl, tmpL):
	''' Function to create string zero string vector before graph filename.
	According to the total number of reactions N zeros will be add before the instant reaction number 
	(e.g. reaction 130 of 10000 the string became '00130')'''
	strZero = ''
	nZeros = len(str(tmpL)) - len(str(tmpl))
	if nZeros > 0:
		for i in range(0,nZeros): strZero = strZero + '0'
	return strZero

if __name__ == '__main__':
	parser = ArgumentParser(
				description='Function to evaluate the activity of each species during the simulation, \
				catalyst substrate product or nothing. Moreover the script recognize all those molecules functioning as hub'
				, epilog='''File with angle trajectories are created. ''') 
	parser.add_argument('-p', '--StrPath', help='Path where files are stored', default='./')
	args = parser.parse_args()
	
	print "Simulation Results Path: ", args.StrPath
	
	StrPath = os.path.abspath(args.StrPath)
		
	tmpDirs = sort(os.listdir(StrPath))
	
	os.chdir(StrPath)
	
	#currentDir = StrPath.split("/")[-1]
	currentDir = ''
	ndn = currentDir + '_0_new_allStatResults'
	newdirAllResults = os.path.abspath(os.path.join(os.curdir, ndn))
	if not os.path.isdir(newdirAllResults):
		try:
			os.mkdir(newdirAllResults)
		except:
			print "Impossible to create statistic directory", newdirAllResults; sys.exit(1)
			
	validDir = 1
	for IDdir, tmpDir in enumerate(tmpDirs):
		dupTime = []
		totDirName = os.path.join(StrPath,tmpDir)
		if os.path.isdir(totDirName):
			# Move to the directory 
			os.chdir(newdirAllResults)
			f_name = "deltat_" + tmpDir +".csv"
			fid_deltat = open(f_name, 'w')
		  	os.chdir(StrPath)			
			os.chdir(totDirName)
			resDirPath = os.path.abspath(os.path.join("./", "res"))
			print " |- Results Folder: ", resDirPath
			if os.path.isdir(resDirPath):
				os.chdir(resDirPath)
				# Find the number of generations
				numberOfGen = len(glob.glob(os.path.join(resDirPath,'times_*')))
								
				for ngen in range(1,numberOfGen+1):
				  
					  strZeros = zeroBeforeStrNum(ngen, numberOfGen)
					  
					  strSpecies = 'timeSpeciesAmount_' + strZeros + str(ngen) + '*'
					  #strSpecies = 'timeSpeciesAmount_00' + strZeros + str(ngen) + '*'  
					  
					   # Searching for files
					  speciesFiles = sorted(glob.glob(os.path.join(resDirPath,strSpecies)))
					  
					  for idS, sngSpeciesFile in enumerate(speciesFiles):
					  	
					  	#print '  |- Species File: ', sngSpeciesFile	
					  	data = np.loadtxt(open(sngSpeciesFile,"rb"),delimiter="\t")
					  	print data[-1,1]
					  	dupTime.append(data[-1,1])
					  	
			  	
			  	
			  	print dupTime
			  	
				# Creare file where store data
				
			  	np.savetxt(fid_deltat, dupTime, fmt='%.4f', delimiter='\t')

				
						
						