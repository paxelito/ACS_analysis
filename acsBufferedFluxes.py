#!/usr/bin/env python
# -*- coding: latin-1 -*-
'''Function to evaluate the activity of each species during the simulation, 
   catalyst substrate product or nothing
'''

import sys, os # Standard library
import datetime as dt
import linecache as lc
import glob
import numpy as np # Scientific library
from numpy import * 
from argparse import ArgumentParser

from lib.IO import *
from lib.graph import network
#--------------------------------------------------------------------------------------
def zeroBeforeStrNum(tmpl, tmpL):
	strZero = ''
	nZeros = len(str(tmpL)) - len(str(tmpl))
	if nZeros > 0:
		for i in range(0,nZeros): strZero = strZero + '0'
	return strZero
#--------------------------------------------------------------------------------------

# Input parameters definition 
if __name__ == '__main__':
	parser = ArgumentParser(
				description='This script characterize the buffered flux dynamics.'
				, epilog='''Buffered flux dynamics. ''') 
	parser.add_argument('-p', '--strPath', help='Path where files are stored', default='./')
	args = parser.parse_args()
	
	# Create absolute paths
	StrPath = os.path.abspath(args.strPath)
	tmpDirs = sort(os.listdir(StrPath))
	
	# Reaction types
	_CONDENSATION_ = 0
	_CLEAVAGE_ = 1
	_ENDOCONDENSATION_ = 7
	_ENDOCLEAVAGE_ = 6
	_SPONTCONDENSATION_ = 10 
	_SPONTCLEAVAGE_ = 11
	
	# Goes into the simulation folder
	os.chdir(StrPath)
	chemistry = 1
	
	# Create stas folders
	ndn = '_0_new_allStatResults'
	newdirAllResults = os.path.join(StrPath, ndn)
	if not os.path.isdir(newdirAllResults):
		try:
			os.mkdir(newdirAllResults)
		except:
			print "Impossible to create statistic directory", newdirAllResults; sys.exit(1)
	
	for tmpDir in tmpDirs:
	
		totDirName = os.path.join(StrPath,tmpDir)
		if os.path.isdir(totDirName):
			# Move to the directory 
			os.chdir(totDirName)
			resDirPath = os.path.abspath(os.path.join("./", "res"))
			print " |- Results Folder: ", resDirPath
			if os.path.isdir(resDirPath):
				os.chdir(resDirPath)
				
				# Find the number of generations
				numberOfGen = len(glob.glob(os.path.join(resDirPath,'times_*')))
				#ÊFor each generation
				os.chdir(resDirPath)
				for ngen in range(1,numberOfGen+1):
					
					tmpFluxFile = 'fluxDynamics_'  + str(tmpDir) + "_" + str(ngen) + '.csv'
					fid_tmpSpeciesStatsSummaryName = open(os.path.join(newdirAllResults,tmpFluxFile), 'w')
					strToWrite = "Time\tIN\tOUT\tDELTA\tbrick_IN\tbrick_OUT\tbrick_DELTA\n"
					fid_tmpSpeciesStatsSummaryName.write(strToWrite)
						
					strZeros = zeroBeforeStrNum(ngen, numberOfGen)
					if ngen == 1:
				  		strSpeciesZero = 'species_' + strZeros + str(0) + '*'
				  		speciesFilesZero = sorted(glob.glob(os.path.join(resDirPath,strSpeciesZero)))
	
					strSpecies = 'species_' + strZeros + str(ngen) + '*'  
					  
				  	# Searching for files
				  	speciesFiles = sorted(glob.glob(os.path.join(resDirPath,strSpecies)))
				  	
				  	# Select the last file to create the specie vector. 
				  	lastfilespecies = speciesFiles[-1]
				  	
				  	print '  |- Create fluxed species vector from: ', lastfilespecies	
					try:
						fidSpecies = open(lastfilespecies, 'r')
					except:
						print ' |- impossible to load ', lastfilespecies; sys.exit(1)
							
					flux_seq = []; nSpecies = 0; flux_ids = []; flux_ssize = [];				
					for sp in fidSpecies:
						tmpID, tmpSeq, tmpConc, tmpDiff, tmpSol, tmpCpxDiss, tmpCpxCut, tmpEval, tmpAge, tmpReb, tmpCatID, tmpSubID, tmpKpho, tmpLoadConc, tmpConcLock = sp.split()
						if int(tmpConcLock) == 1:
							flux_seq.append(str(tmpSeq))
							flux_ids.append(int(tmpID))
							flux_ssize.append(int(len(str(tmpSeq))))
							nSpecies += 1
						
					# Create lists containing stats data 
					# The number of columns is equal to times + tot(MOLS) + tot(Bricks) + each buffered species (mols and bricks)
					# counters = np.zeros((nSpecies,(len(flux_seq)*2)+6)) 
					# Select reaction_parameter file
					strRctPar = 'reactions_parameters_' + strZeros + str(ngen) + '_1.csv'
					
				  	print '  |- Open reactions parameter file: : ', strRctPar	
					try:
						fidRctPar = open(strRctPar, 'r')
					except:
						print ' |- impossible to load ', strRctPar; sys.exit(1)	
						
					totIN = 0
					totOUT = 0
					deltaIO = 0
					tot_B_IN = 0
					tot_B_OUT = 0
					delta_B_IO = 0
					
					strRct = 'reactions_' + strZeros + str(ngen) + '*'  
					strCat = 'catalysis_' + strZeros + str(ngen) + '*'  
				  	# Searching for files
				  	rctFiles = sorted(glob.glob(os.path.join(resDirPath,strRct)))
				  	catFiles = sorted(glob.glob(os.path.join(resDirPath,strCat)))
				  	# Upload latest file
				  	lastRct = readfiles.loadAllData(totDirName,rctFiles[-1]) # reaction file upload
				  	lastCat = readfiles.loadAllData(totDirName,catFiles[-1]) # reaction file upload
					
					for idRct, sngRct in enumerate(fidRctPar): 
						# Load single reaction parameters
						
						if idRct % 10000 == 0: print "\t|- Reaction: ", idRct
						
						reaction, rctTime, rctType, cat, S1, S2, S3, loadedMolsConc, loadedMols,\
   						gillMean, gillSD, gillEntropy, newSpeciesCreationProb, reverseProbability = readfiles.splitRctParsLine(sngRct)
   						
   						#S1, S2, S3 = network.fixCondensationReaction(S1, S2, S3, lastRct)
						
						# counters[idRct,0] = rctTime
						
						if (rctType == _CONDENSATION_) | (rctType == _ENDOCONDENSATION_) | (rctType == _SPONTCONDENSATION_):
							if (S1 in flux_ids): 
								deltaIO -= 1
								totOUT += 1
								tot_B_OUT += len(flux_seq[S1])
								delta_B_IO -= len(flux_seq[S1])
							if (S2 in flux_ids): 
								totIN += 1
								deltaIO += 1
								tot_B_IN += len(flux_seq[S2])
								delta_B_IO += len(flux_seq[S2])
							if (S3 in flux_ids): 
								totIN += 1	
								deltaIO += 1
								tot_B_IN += len(flux_seq[S3])
								delta_B_IO += len(flux_seq[S3])
						else:
							if (S1 in flux_ids): 
								deltaIO += 1
								totIN += 1
								tot_B_IN += len(flux_seq[S1])
								delta_B_IO += len(flux_seq[S1])
							if (S2 in flux_ids): 
								totOUT += 1
								deltaIO -= 1
								tot_B_OUT += len(flux_seq[S2])
								delta_B_IO -= len(flux_seq[S2])
							if (S3 in flux_ids): 
								totOUT += 1	
								deltaIO -= 1
								tot_B_OUT += len(flux_seq[S3])
								delta_B_IO -= len(flux_seq[S3])
								
						str2wrt = str(rctTime) + "\t" + str(totIN) + "\t" + str(totOUT) + "\t" + str(deltaIO) + "\t" + str(tot_B_IN) + "\t" + str(tot_B_OUT) + "\t"  + str(delta_B_IO) + "\n"
						fid_tmpSpeciesStatsSummaryName.write(str2wrt) 
						
					fid_tmpSpeciesStatsSummaryName.close()
					
	
					
					
			chemistry += 1	
					
				  	
				  	


