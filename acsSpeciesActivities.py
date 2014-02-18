#!/usr/bin/env python
# -*- coding: latin-1 -*-
'''Function to evaluate the activity of each species during the simulation, 
   catalyst substrate product or nothing. Moreover the script recognize all those molecules functioning as hub
'''

import sys, os # Standard library
import datetime as dt
import linecache as lc
import glob
import numpy as np # Scientific library
from numpy import * 

try:
    from pylab import *
except:
    pass
   
#--------------------------------------------------------------------------------------
def zeroBeforeStrNum(tmpl, tmpL):
	strZero = ''
	nZeros = len(str(tmpL)) - len(str(tmpl))
	if nZeros > 0:
		for i in range(0,nZeros): strZero = strZero + '0'
	return strZero
#--------------------------------------------------------------------------------------
	
if __name__ == '__main__':
	parser = ArgumentParser(
				description='Function to evaluate the activity of each species during the simulation, \
				catalyst substrate product or nothing. Moreover the script recognize all those molecules functioning as hub'
				, epilog='''File with angle trajectories are created. ''') 
	parser.add_argument('-p', '--strPath', help='Path where files are stored', default='./')
	args = parser.parse_args()


	StrPath = os.path.abspath(args.StrPath)
		
	tmpDirs = sort(os.listdir(StrPath))
	
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
	
	ndn = 'allSpeciesActivities'
	newdirAllResultsInner = os.path.join(StrPath,'_0_new_allStatResults',ndn)
	if not os.path.isdir(newdirAllResultsInner):
		try:
			os.mkdir(newdirAllResultsInner)
		except:
			print "Impossible to create all species acts directory", newdirAllResults; sys.exit(1)
	
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
				# Go into stat folders
				os.chdir(newdirAllResults)	
				os.chdir(newdirAllResultsInner)	
				tmpSpeciesStatsSummaryName = 'speciesStatsSummary_'  + str(chemistry) + '.csv'
				tmpSpeciesStatsSummaryNameFID = open(tmpSpeciesStatsSummaryName, 'w')
				# Change Chemistry
				tmpStr = '\n------- CHEMISTRY ' + str(chemistry) + '---------\n'
				tmpSpeciesStatsSummaryNameFID.write(tmpStr)
				os.chdir(resDirPath)
				for ngen in range(1,numberOfGen+1):
					
					tmpStr = '\n- - - - Generation ' + str(ngen) + '- - - - -\n'
					tmpSpeciesStatsSummaryNameFID.write(tmpStr)
						
					strZeros = zeroBeforeStrNum(ngen, numberOfGen)
					if ngen == 1:
				  		strSpeciesZero = 'species_' + strZeros + str(0) + '*'
				  		speciesFilesZero = sorted(glob.glob(os.path.join(resDirPath,strSpeciesZero)))
	
					strSpecies = 'species_' + strZeros + str(ngen) + '*'  
					  
				  	# Searching for files
				  	speciesFiles = sorted(glob.glob(os.path.join(resDirPath,strSpecies)))
				  	
				  	# Select the last file to create the specie vector. 
				  	lastfilespecies = speciesFiles[-1]
				  	
				  	print '  |- Create species vector from: ', lastfilespecies	
					try:
						fidSpecies = open(lastfilespecies, 'r')
					except:
						print ' |- impossible to load ', lastfilespecies; sys.exit(1)
							
					seq = []; nSpecies = 0
					for sp in fidSpecies:
						tmpID, tmpSeq, tmpConc, tmpDiff, tmpSol, tmpCpxDiss, tmpCpxCut, tmpEval, tmpAge, tmpReb, tmpCatID, tmpSubID, tmpKpho, tmpLoadConc, tmpConcLock = sp.split()
						seq.append(str(tmpSeq))
						nSpecies += 1
						
					# Create lists containing statistical data 
					counters = np.zeros((nSpecies,3)) 
					# Select reaction_parameter file
					strRctPar = 'reactions_parameters_' + strZeros + str(ngen) + '_1.csv'
					
				  	print '  |- Open reactions parameter file: : ', strRctPar	
					try:
						fidRctPar = open(strRctPar, 'r')
					except:
						print ' |- impossible to load ', strRctPar; sys.exit(1)	
					
					for idRct, sngRct in enumerate(fidRctPar): 
						tmpRctId, tmpRctT, tmpRctType, tmpCat, tmpS1, tmpS2, tmpS3, tmpCM, \
						tmpCC, tmpGil, tmpGilSD, tmpE, tmpNSP, tmpRatio = sngRct.split()
						rctType = int(tmpRctType)
						cat = int(tmpCat)
						S1 = int(tmpS1)
						S2 = int(tmpS2)
						S3 = int(tmpS3)
						if (rctType == 0) | (rctType == 7):
							#catalyst
							counters[cat,0] += 1
							#sub1
							counters[S2,1] += 1
							#sub2
							counters[S3,1] += 1
							#prod
							counters[S1,2] += 1
						else:
							#catalyst
							counters[cat,0] += 1
							#sub1
							counters[S2,1] += 1
							#prod1
							counters[S3,2] += 1
							#prod2
							counters[S1,2] += 1		
					
					
					# Go into stat folders
					os.chdir(newdirAllResults)	
					os.chdir(newdirAllResultsInner)		
					
					
					tmpFileName = 'speciesStats_'  + str(chemistry) + '_' + strZeros + str(ngen) + '.csv'
					print '  |- Saving result file: : ', tmpFileName
					tmpFileNameFID = open(tmpFileName, 'w')
					ID = 0
					tmpStr = 'Total Number of Reactions\t\t\t\t' + str(idRct) + '\n'
					tmpSpeciesStatsSummaryNameFID.write(tmpStr)
					tmpStr = 'ID\tCat\tSub\tProd\n'
					tmpSpeciesStatsSummaryNameFID.write(tmpStr)
					for sngCnt in counters:
						tmpStr = str(ID) + "\t" + str(sngCnt[0]) + "\t" + str(sngCnt[1]) + "\t" + str(sngCnt[2]) + "\n"
						tmpFileNameFID.write(tmpStr)
						if (sngCnt[0] > (sngCnt[1] * 100)) & (sngCnt[0] > 1000):
							tmpSpeciesStatsSummaryNameFID.write(tmpStr)
						ID += 1
					tmpFileNameFID.close()
					
					os.chdir(resDirPath)	
				
				tmpSpeciesStatsSummaryNameFID.close()
				chemistry += 1	
					
				  	
			  	


