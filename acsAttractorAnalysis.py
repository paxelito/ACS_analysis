#!/usr/bin/env python
'''Function to compute the different attractors emerging from a simulazione. The algorithm compares all the final states of the simulation computing 
	The differences between those
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
	
# get path for placing simulation
try:
	StrPath = sys.argv[1] # Here the path of the simulation output file
	tmpMaxFluxL = int(sys.argv[2]) # Influx max length
except:
	print "Usage:",sys.argv[0], "infile outfile"; sys.exit(1)
	
print "Simulation Results Path: ", StrPath

today = dt.date.today()
	
tmpDirs = os.listdir(StrPath)
allSortedSpecies = [] 
allConcentrations = [] 
allSortedSpeciesNOINFLUX = [] 
allConcentrationsNOINFLUX = [] 
print ''
print 'ACS ATTRACTORS ANALYSER'
print ''
print '|- STEP 1. Creating common sorted species list...'
for tmpDir in tmpDirs:
	os.chdir(StrPath)
	totDirName = os.path.join(StrPath,tmpDir)
	if os.path.isdir(totDirName):
		# Move to the directory 
		os.chdir(totDirName)
		resDirPath = os.path.join(totDirName, "res")
		if os.path.isdir(resDirPath):
			os.chdir(resDirPath)
		
			# Searching for files
			speciesFiles = sorted(glob.glob(os.path.join(resDirPath,'species_*')))
			speciesFile = speciesFiles[-1]
		
			print " |- Species File: ", speciesFile
			# Open Catalysis File
			try:
				fidSpecies = open(speciesFile, 'r')
			except:
				print ' |- impossible to load ', speciesFile; sys.exit(1)
			
			# For each last species file
			seq = []; conc = []; seqNOINFLUX = []
			for sp in fidSpecies:
				tmpID, tmpSeq, tmpConc, tmpDiff, tmpSol, tmpCpxDiss, tmpCpxCut, tmpEval, tmpAge, tmpReb, tmpCatID, tmpSubID, tmpKpho, tmpLoadConc, tmpConcLock = sp.split()
				if (int(tmpCpxCut) == 0) & (float(tmpConc) > 0):
					seq.append(str(tmpSeq))
					conc.append(float(tmpConc))
					if len(str(tmpSeq)) > tmpMaxFluxL:
						seqNOINFLUX.append(str(tmpSeq))	
			# Common ordered species List 
			allSortedSpecies = list(set(seq) | set(allSortedSpecies))
			allSortedSpecies.sort()	
			# NO influx
			allSortedSpeciesNOINFLUX = list(set(seqNOINFLUX) | set(allSortedSpeciesNOINFLUX))
			allSortedSpeciesNOINFLUX.sort()	
			
			fidSpecies.close()	
			
print '|- STEP 2. Creating species concentration lists according to the overall sorted species list...'		
overallConcList = []
overallConcListNOINFLUX = []
numberOfFolders = 0
for tmpDir in tmpDirs:
	os.chdir(StrPath)
	totDirName = os.path.join(StrPath,tmpDir)
	if os.path.isdir(totDirName):
		# Move to the directory 
		os.chdir(totDirName)
		resDirPath = os.path.join(totDirName, "res")
		if os.path.isdir(resDirPath):
			numberOfFolders += 1
			os.chdir(resDirPath)
		
			# Searching for files
			speciesFiles = sorted(glob.glob(os.path.join(resDirPath,'species_*')))
			speciesFile = speciesFiles[-1]
		
			print " |- Species File: ", speciesFile
			# Open Catalysis File
			try:
				fidSpecies = open(speciesFile, 'r')
			except:
				print ' |- impossible to load ', speciesFile; sys.exit(1)
			
			# For each last species file
			seq = []; conc = []; speciesConc = []; speciesConcNOINFLUX = []
			for sp in fidSpecies:
				tmpID, tmpSeq, tmpConc, tmpDiff, tmpSol, tmpCpxDiss, tmpCpxCut, tmpEval, tmpAge, tmpReb, tmpCatID, tmpSubID, tmpKpho, tmpLoadConc, tmpConcLock = sp.split()
				if (int(tmpCpxCut) == 0) & (float(tmpConc) > 0):
					seq.append(str(tmpSeq))
					conc.append(float(tmpConc))
					
			# Check the presence of the species in the all common species list		
			for key in allSortedSpecies:
				try:
					pos = seq.index(key)
					speciesConc.append(conc[pos])
				except:
					speciesConc.append(0)
			for key2 in allSortedSpeciesNOINFLUX:
				try:
					pos2 = seq.index(key2)
					speciesConcNOINFLUX.append(conc[pos2])
				except:
					speciesConcNOINFLUX.append(0)
			
			# Add the concentration list in the concentrationS list
			overallConcList.append(speciesConc)
			overallConcListNOINFLUX.append(speciesConcNOINFLUX)		
			
			fidSpecies.close()	
			
print '|- STEP 3. Compute attractors differences (in term of different multi-dimensional angles)'	
overallResMatrix = np.zeros((numberOfFolders,numberOfFolders))
overallResMatrix_angle = np.zeros((numberOfFolders,numberOfFolders))
overallResMatrixNOINFLUX = np.zeros((numberOfFolders,numberOfFolders))
overallResMatrix_angleNOINFLUX = np.zeros((numberOfFolders,numberOfFolders))
# Angle between lists is now computed
for idx, lx in enumerate(overallConcList):
	for idy, ly in enumerate(overallConcList):
		vecX = np.array(lx)
		vecY = np.array(ly)
		# Compute coseno
		tmpCos = np.dot(vecX,vecY) / (np.linalg.norm(vecX) * np.linalg.norm(vecY))
		overallResMatrix[idx,idy] = tmpCos
		overallResMatrix_angle[idx,idy] = np.arccos(tmpCos)
for idx, lx in enumerate(overallConcListNOINFLUX):
	for idy, ly in enumerate(overallConcListNOINFLUX):
		vecX = np.array(lx)
		vecY = np.array(ly)
		# Compute coseno
		tmpCos = np.dot(vecX,vecY) / (np.linalg.norm(vecX) * np.linalg.norm(vecY))
		overallResMatrixNOINFLUX[idx,idy] = tmpCos
		overallResMatrix_angleNOINFLUX[idx,idy] = np.arccos(tmpCos)

print '|- STEP 3. Save Files'
os.chdir(StrPath)
# '''Function to save statistic on file'''
outFnameStat = 'acsAttractorsAnalysis.csv'
saveFileStat = open(outFnameStat, 'w')
cnt = 0
for i in range(numberOfFolders):
	strTypes = ''
	for j in range(numberOfFolders):
		strTypes += str(overallResMatrix[j,i]) + '\t'	
	strTypes += '\n'
	saveFileStat.write(strTypes)
	cnt += 1
saveFileStat.close()	

outFnameStat = 'acsAttractorsAnalysisNOINFLUX.csv'
saveFileStat = open(outFnameStat, 'w')
cnt = 0
for i in range(numberOfFolders):
	strTypes = ''
	for j in range(numberOfFolders):
		strTypes += str(overallResMatrixNOINFLUX[j,i]) + '\t'	
	strTypes += '\n'
	saveFileStat.write(strTypes)
	cnt += 1
saveFileStat.close()	



print '|- FINISHED... SEE YOU NEXT TIME'		
			

			
					
						
				
			
			

			
    			


					
					

			
				

	
	
