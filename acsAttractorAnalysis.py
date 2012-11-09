#!/usr/bin/env python
'''Function to compute the different attractors emerging from a simulation. The algorithm compares all the final states of the simulation computing 
	The differences between those. test
	https://help.github.com/articles/fork-a-repo
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
# Function to create string zero string vector before graph filename. 
# According to the total number of reactions N zeros will be add before the instant reaction number 
# (e.g. reaction 130 of 10000 the string became '00130')
def zeroBeforeStrNum(tmpl, tmpL):
	strZero = ''
	nZeros = len(str(tmpL)) - len(str(tmpl))
	if nZeros > 0:
		for i in range(0,nZeros): strZero = strZero + '0'
	return strZero
	
# get path for placing simulation
try:
	StrPath = sys.argv[1] # Here the path of the simulation output file
	tmpMaxFluxL = int(sys.argv[2]) # Influx max length
except:
	print "Usage:",sys.argv[0], "infile outfile"; sys.exit(1)
	
print "Simulation Results Path: ", StrPath

# Convert path in absolute path
StrPath = os.path.abspath(StrPath)

today = dt.date.today()

StrPath = os.path.abspath(StrPath)
	
tmpDirs = sort(os.listdir(StrPath))
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
		resDirPath = os.path.abspath(os.path.join("./", "res"))
		print " |- Results Folder: ", resDirPath
		if os.path.isdir(resDirPath):
			os.chdir(resDirPath)
			
			# Find the number of generations
			numberOfGen = len(glob.glob(os.path.join(resDirPath,'times_*')))
			
			for ngen in range(1,numberOfGen+1):
			  
			  strZeros = zeroBeforeStrNum(ngen, numberOfGen)
			  strSpecies = 'species_' + strZeros + str(ngen) + '*'
				  
			  # Searching for files
			  speciesFiles = sorted(glob.glob(os.path.join(resDirPath,strSpecies)))
			  speciesFile = speciesFiles[-1]
		  
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
		else: 
			  print " |- no result folder has been found"
			
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
			
			os.chdir(resDirPath)
			
			# Find the number of generations
			numberOfGen = len(glob.glob(os.path.join(resDirPath,'times_*')))
		
			for ngen in range(1,numberOfGen+1):
			  
			  numberOfFolders += 1
			  
			  strZeros = zeroBeforeStrNum(ngen, numberOfGen)
			  strSpecies = 'species_' + strZeros + str(ngen) + '*'
		  
			  # Searching for files
			  speciesFiles = sorted(glob.glob(os.path.join(resDirPath,strSpecies)))
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
		#overallResMatrix_angle[idx,idy] = np.arccos(tmpCos)
for idx, lx in enumerate(overallConcListNOINFLUX):
	for idy, ly in enumerate(overallConcListNOINFLUX):
		vecX = np.array(lx)
		vecY = np.array(ly)
		# Compute coseno
		tmpCos = np.dot(vecX,vecY) / (np.linalg.norm(vecX) * np.linalg.norm(vecY))
		overallResMatrixNOINFLUX[idx,idy] = tmpCos
		#overallResMatrix_angleNOINFLUX[idx,idy] = np.arccos(tmpCos)

print '|- STEP 3. Save Files'
os.chdir(StrPath)
# ''' Save concentration vectors"
print ' |- Saving species conentrations list...'
outFnameStat = 'speciesVector.csv'
saveFileStat = open(outFnameStat, 'w')

for sngSimConc in overallConcList:
	strTypes = ''
	for j in allSortedSpecies:
		strTypes += j + '\t'
	strTypes += '\n'
	for j in sngSimConc:
		strTypes += str(j) + '\t'
	strTypes += '\n'
	saveFileStat.write(strTypes)
saveFileStat.close()

# ''' Save concentration vectors (NO INFLUX)"
print ' |- Saving species conentrations list (NO INFLUX)...'
outFnameStat = 'speciesVectorNOINFLUX.csv'
saveFileStat = open(outFnameStat, 'w')
for id, sngSimConc in enumerate(overallConcListNOINFLUX):
	strTypes = ''
	for j in sngSimConc:
		strTypes += allSortedSpeciesNOINFLUX[id] + '\t' + str(j) + '\t'
	strTypes += '\n'
	saveFileStat.write(strTypes)
saveFileStat.close()

# '''Function to save statistic on file'''
print ' |- Saving coseno...'
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

print ' |- Saving coseno (NO INFLUX)...'
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
			

			
					
						
				
			
			

			
    			


					
					

			
				

	
	
