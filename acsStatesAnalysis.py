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

def angleBetweenTwoLists(tmpSeqX, tmpConcX, tmpSeqY, tmpConcY):
	tmpAllList = list(set(tmpSeqX) | set(tmpSeqY))
	tmpAllList.sort()	
	speciesConcX = []
	speciesConcY = []
	for key in tmpAllList:
		try:
			posX = tmpSeqX.index(key)
			speciesConcX.append(tmpConcX[posX])
		except:
			speciesConcX.append(0)
		try:
			posY = tmpSeqY.index(key)
			speciesConcY.append(tmpConcY[posY])
		except:
			speciesConcY.append(0)
	 
	vecX = np.array(speciesConcX)
	vecY = np.array(speciesConcY)
	# Compute coseno AND WRITE FILE
	tmpCos = np.dot(vecX,vecY) / (np.linalg.norm(vecX) * np.linalg.norm(vecY))
	
	return tmpCos
	
# get path for placing simulation
try:
	StrPath = sys.argv[1] # Here the path of the simulation output file
	tmpMaxFluxL = int(sys.argv[2]) # Influx max length
except:
	print "Usage:",sys.argv[0], "infile outfile"; sys.exit(1)
	
print "Simulation Results Path: ", StrPath

today = dt.date.today()

StrPath = os.path.abspath(StrPath)
	
tmpDirs = sort(os.listdir(StrPath))

os.chdir(StrPath)
 
print ''
print 'ACS STATES ANALYSER'
print ''
print '|- STEP 1. Creating common sorted species list...'
# Open File containing results 
previousFILEname = 't_tminus_1.csv'
previousFILE_FID = open(previousFILEname, 'w')
previousNOINFLUX_FILEname = 't_tminus_1_NOINFLUX.csv'
previousNOINFLUX_FILE_FID = open(previousNOINFLUX_FILEname, 'w')
startFILEname = 't_start.csv'
startFILE_FID = open(startFILEname, 'w')
startNOINFLUX_FILEname = 't_start_NOINFLUX.csv'
startNOINFLUX_FILE_FID = open(startNOINFLUX_FILEname, 'w')
newSpeciesFileName = 'newSpecies.csv'
newSpecies_FID = open(newSpeciesFileName, 'w')
newSpeciesFileName = 'livingSpecies.csv'
livingSpecies_FID = open(newSpeciesFileName, 'w')
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
			
			for ngen in range(1,numberOfGen+1):
			  
				  strZeros = zeroBeforeStrNum(ngen, numberOfGen)
				  
				  if ngen == 1:
				  	strSpeciesZero = 'species_' + strZeros + str(0) + '*'
				  	speciesFilesZero = sorted(glob.glob(os.path.join(resDirPath,strSpeciesZero)))
				  				  
				  strSpecies = 'species_' + strZeros + str(ngen) + '*'  
					  
				  # Searching for files
				  speciesFiles = sorted(glob.glob(os.path.join(resDirPath,strSpecies)))
				  
				  if ngen == 1:
				  	speciesFiles = speciesFilesZero + speciesFiles
				  	
				  # FOR EACH FILE SPECIES
				  seqOLD = []; seqOLDNOINFLUX = []; concOLD = []
				  previousAngleList = []
				  previousAngleListNOINFLUX = []
				  seqSTART = []; seqSTART_NOINFLUX = []; concSTART = []
				  oldNumberOfSpecies = 0
				  for idS, sngSpeciesFile in enumerate(speciesFiles):
				  	
				  	print '  |- Species File: ', sngSpeciesFile
						
					# Open Catalysis File
					try:
						fidSpecies = open(sngSpeciesFile, 'r')
					except:
						print ' |- impossible to load ', sngSpeciesFile; sys.exit(1)
						
					# For each last species file
					seq = []; conc = []; seqNOINFLUX = []; numberOfSpecies = 0;
					for sp in fidSpecies:
						tmpID, tmpSeq, tmpConc, tmpDiff, tmpSol, tmpCpxDiss, tmpCpxCut, tmpEval, tmpAge, tmpReb, tmpCatID, tmpSubID, tmpKpho, tmpLoadConc, tmpConcLock = sp.split()
						if (int(tmpCpxCut) == 0) & (float(tmpConc) > 0):
							seq.append(str(tmpSeq))
							conc.append(float(tmpConc))
							if idS == 0:
								seqSTART.append(str(tmpSeq))
								concSTART.append(float(tmpConc))							
							if len(str(tmpSeq)) > tmpMaxFluxL:
								seqNOINFLUX.append(str(tmpSeq))	
								if idS == 0:
									seqSTART_NOINFLUX.append(str(tmpSeq))
						if int(tmpCpxCut) == 0:
							numberOfSpecies += 1	
							
					deltaNspecies = numberOfSpecies - oldNumberOfSpecies
					oldNumberOfSpecies = numberOfSpecies 
					strtoW = str(deltaNspecies) + '\t'
					newSpecies_FID.write(strtoW)	
					strtoW = str(len(conc)) + '\t'
					livingSpecies_FID.write(strtoW)						
					
					# ------------------------------------------------------					
					# PREVIOUS ONE Defining concentration of the two vectors
					coseno = angleBetweenTwoLists(seq, conc, seqOLD, concOLD)
					if idS != 0:
						previousAngleList.append(coseno)
						strtoW = str(coseno) + '\t'
						previousFILE_FID.write(strtoW)

					# ------------------------------------------------------					
					# START Defining concentration of the two vectors
					coseno = angleBetweenTwoLists(seq, conc, seqSTART, concSTART)
					if idS != 0:
						previousAngleList.append(coseno)
						strtoW = str(coseno) + '\t'
						startFILE_FID.write(strtoW)						
						  
					# ------------------------------------------------------					
					# PREVIOUS ONE (NO INFLUX) Defining concentration of the two vectors
					coseno = angleBetweenTwoLists(seqNOINFLUX, conc, seqOLDNOINFLUX, concOLD)
					if idS != 0:
						previousAngleListNOINFLUX.append(coseno)	
						strtoW = str(coseno) + '\t'
						previousNOINFLUX_FILE_FID.write(strtoW)		

					# ------------------------------------------------------					
					# START (NO INFLUX) Defining concentration of the two vectors
					coseno = angleBetweenTwoLists(seqNOINFLUX, conc, seqSTART_NOINFLUX, concSTART)
					if idS != 0:
						previousAngleList.append(coseno)
						strtoW = str(coseno) + '\t'
						startNOINFLUX_FILE_FID.write(strtoW)								
						  
					# the new lists becomes the old one
					seqOLD = seq[:]
					seqOLDNOINFLUX = seqNOINFLUX[:]
					concOLD = conc[:]
							  					  
				  	fidSpecies.close()
				  	
				  previousFILE_FID.write('\n')
				  previousNOINFLUX_FILE_FID.write('\n')
				  startFILE_FID.write('\n')
				  startNOINFLUX_FILE_FID.write('\n')
				  newSpecies_FID.write('\n')
				  livingSpecies_FID.write('\n')
		else: 
			print " |- no result folder has been found"

# CLOSE FILES
previousFILE_FID.close()
previousNOINFLUX_FILE_FID.close()
startFILE_FID.close()
startNOINFLUX_FILE_FID.close()
newSpecies_FID.close()
livingSpecies_FID.close()


print '|- FINISHED... SEE YOU NEXT TIME'		
			

			
					
						
				
			
			

			
    			


					
					

			
				

	
	
