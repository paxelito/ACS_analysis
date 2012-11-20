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

def zeroBeforeStrNum(tmpl, tmpL):
	''' Function to create string zero string vector before graph filename.
	According to the total number of reactions N zeros will be add before the instant reaction number 
	(e.g. reaction 130 of 10000 the string became '00130')'''
	strZero = ''
	nZeros = len(str(tmpL)) - len(str(tmpl))
	if nZeros > 0:
		for i in range(0,nZeros): strZero = strZero + '0'
	return strZero

def returnZeroSpeciesList(tmpLastSpeciesFile):
	'''Function to create a zero vector for each species (NO COMPLEXES)'''
	try:
		fidSpecies = open(tmpLastSpeciesFile, 'r')
	except:
		print ' |- impossible to load ', tmpLastSpeciesFile; sys.exit(1)
	tmpZeroList = []
	for s in fidSpecies:
		tmpID, tmpSeq, tmpConc, tmpDiff, tmpSol, tmpCpxDiss, tmpCpxCut, tmpEval, tmpAge, tmpReb, tmpCatID, tmpSubID, tmpKpho, tmpLoadConc, tmpConcLock = s.split()
		if (int(tmpCpxCut) == 0):
			tmpZeroList.append(0)
		
	return tmpZeroList

def distanceMisures(tmpSeqX, tmpConcX, tmpSeqY, tmpConcY, tmpIDs):
	'''Function to compute the angle between two multidimensional vectors'''
	strtoW = [0,0,0]
	if tmpIDs != 0:
		# COS BETWEEN VECTORS
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
		tmpCos = np.dot(vecX,vecY) / (np.linalg.norm(vecX) * np.linalg.norm(vecY))	
		
		strtoW[0] = tmpCos
		
		# HAMMING DISTANCE and EUCLIDEAN DISTANCE
		tmpHD = 0
		tmpEU = 0
		for pos, x in enumerate(speciesConcX):
			if ((x > 0.0) & (speciesConcY[pos] == 0.0)) | ((x == 0.0) & (speciesConcY[pos] > 0.0)):
				tmpHD += 1
			tmpEU += pow(x - speciesConcY[pos],2)
			
		strtoW[1] = tmpHD
		strtoW[2] = pow(tmpEU,0.5)
					
	else:
		strtoW[0] = 1
		strtoW[1] = 0
		strtoW[2] = 0
		
	return strtoW
	
	
try:
	StrPath = sys.argv[1] # Here the path of the simulation output file
	tmpMaxFluxL = int(sys.argv[2]) # Influx max length
	tmpVolume = float(sys.argv[3])  
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

ndn = '_0_new_allStatResults_'
newdirAllResults = os.path.join(os.curdir, ndn)
if not os.path.isdir(newdirAllResults):
	try:
		os.mkdir(newdirAllResults)
	except:
		print "Impossible to create statistic directory", newdirAllResults; sys.exit(1)
		
os.chdir(newdirAllResults)

# Open File containing results 
previousFILE_FID = open('STAT_t_tminus_1.csv', 'w')
previousNOINFLUX_FILE_FID = open('STAT_t_tminus_1_NOINFLUX.csv', 'w')
startFILE_FID = open('STAT_t_start.csv', 'w')
startNOINFLUX_FILE_FID = open('STAT_t_start_NOINFLUX.csv', 'w')

HAM_previousFILE_FID = open('STAT_HAM_t_tminus_1.csv', 'w')
HAM_previousNOINFLUX_FILE_FID = open('STAT_HAM_t_tminus_1_NOINFLUX.csv', 'w')
HAM_startFILE_FID = open('STAT_HAM_t_start.csv', 'w')
HAM_startNOINFLUX_FILE_FID = open('STAT_HAM_t_start_NOINFLUX.csv', 'w')

EUC_previousFILE_FID = open('STAT_EUC_t_tminus_1.csv', 'w')
EUC_previousNOINFLUX_FILE_FID = open('STAT_EUC_t_tminus_1_NOINFLUX.csv', 'w')
EUC_startFILE_FID = open('STAT_EUC_t_start.csv', 'w')
EUC_startNOINFLUX_FILE_FID = open('STAT_EUC_t_start_NOINFLUX.csv', 'w')

previousFILE_FID_group = open('STAT_t_tminus_1_group.csv', 'w')
previousNOINFLUX_FILE_FID_group = open('STAT_t_tminus_1_NOINFLUX_group.csv', 'w')
startFILE_FID_group = open('STAT_t_start_group.csv', 'w')
startNOINFLUX_FILE_FID_group = open('STAT_t_start_NOINFLUX_group.csv', 'w')

HAM_previousFILE_FID_group = open('STAT_HAM_t_tminus_1_group.csv', 'w')
HAM_previousNOINFLUX_FILE_FID_group = open('STAT_HAM_t_tminus_1_NOINFLUX_group.csv', 'w')
HAM_startFILE_FID_group = open('STAT_HAM_t_start_group.csv', 'w')
HAM_startNOINFLUX_FILE_FID_group = open('STAT_HAM_t_start_NOINFLUX_group.csv', 'w')

EUC_previousFILE_FID_group = open('STAT_EUC_t_tminus_1_group.csv', 'w')
EUC_previousNOINFLUX_FILE_FID_group = open('STAT_EUC_t_tminus_1_NOINFLUX_group.csv', 'w')
EUC_startFILE_FID_group = open('STAT_EUC_t_start_group.csv', 'w')
EUC_startNOINFLUX_FILE_FID_group = open('STAT_EUC_t_start_NOINFLUX_group.csv', 'w')

newSpecies_FID = open('STAT_newSpecies.csv', 'w')
livingSpecies_FID = open('STAT_livingSpecies.csv', 'w')
totMass_FID = open('STAT_overallMass.csv', 'w')
evaluatedFID = open('STAT_evaluated.csv', 'w')
zeroOneSpeciesFID = open('STAT_zeroOneSpecies.csv', 'w')

os.chdir(StrPath)

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
			
			group_A_prev = []; group_HAM_prev = []; group_EUC_prev = [];
			group_A_start = []; group_HAM_start = []; group_EUC_start = [];
			group_A_prev_NI = []; group_HAM_prev_NI = []; group_EUC_prev_NI = [];
			group_A_start_NI = []; group_HAM_start_NI = []; group_EUC_start_NI = [];
			
			for ngen in range(1,numberOfGen+1):
			  
				  strZeros = zeroBeforeStrNum(ngen, numberOfGen)
				  
				  if ngen == 1:
				  	strSpeciesZero = 'species_' + strZeros + str(0) + '*'
				  	speciesFilesZero = sorted(glob.glob(os.path.join(resDirPath,strSpeciesZero)))
				  				  
				  strSpecies = 'species_' + strZeros + str(ngen) + '*'  
					  
				  # Searching for files
				  speciesFiles = sorted(glob.glob(os.path.join(resDirPath,strSpecies)))
				  
				  zeroList = returnZeroSpeciesList(speciesFiles[-1])
				  
				  if ngen == 1:
				  	speciesFiles = speciesFilesZero + speciesFiles
				  	
				  # Initialize moving average lists
				  seqOLD = []; seqOLDNOINFLUX = []; concOLD = []
				  seqSTART = []; seqSTART_NOINFLUX = []; concSTART = []
				  totMass = []; obsSpecies = [];
				  
				  oldNumberOfSpecies = 0
				  
				  # FOR EACH FILE SPECIES ---------------------------
				  for idS, sngSpeciesFile in enumerate(speciesFiles):
				  	
				  	print '  |- Species File: ', sngSpeciesFile	
					try:
						fidSpecies = open(sngSpeciesFile, 'r')
					except:
						print ' |- impossible to load ', sngSpeciesFile; sys.exit(1)
						
					seq = []; conc = []; seqNOINFLUX = []; numberOfSpecies = 0; tmpMass = 0; tmpObsSpecies = 0
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
							# Set 1 in zeroList (species has been created at least once)
							zeroList[numberOfSpecies] = 1
							# Update systems mass
							tmpMass += len(str(tmpSeq)) * int(round(float(tmpConc) * 6.022e23 * tmpVolume))
						# If the species is not a complex, the number of species is updated	
						if int(tmpCpxCut) == 0:
							numberOfSpecies += 1	
						if int(tmpCpxCut) == 0 and int(tmpEval) == 1:
							tmpObsSpecies += 1
					
					# Compute number of new species		
					deltaNspecies = numberOfSpecies - oldNumberOfSpecies
					oldNumberOfSpecies = numberOfSpecies 
					strtoW = str(deltaNspecies) + '\t'
					newSpecies_FID.write(strtoW)	
					strtoW = str(len(conc)) + '\t'
					livingSpecies_FID.write(strtoW)		
					strtoW = str(tmpMass) + '\t'
					totMass_FID.write(strtoW)	
					strtoW = str(tmpObsSpecies) + '\t'
					evaluatedFID.write(strtoW)	
					
					# ------------------------------------------------------					
					# PREVIOUS ONE Defining concentration of the two vectors
					tmpMisure = distanceMisures(seq, conc, seqOLD, concOLD, idS)
					previousFILE_FID.write(str(tmpMisure[0]) + '\t'); group_A_prev.append(tmpMisure[0])
					HAM_previousFILE_FID.write(str(tmpMisure[1]) + '\t'); group_HAM_prev.append(tmpMisure[1])
					EUC_previousFILE_FID.write(str(tmpMisure[2]) + '\t'); group_EUC_prev.append(tmpMisure[2])			
					# START Defining concentration of the two vectors
					tmpMisure = distanceMisures(seq, conc, seqSTART, concSTART, idS)
					startFILE_FID.write(str(tmpMisure[0]) + '\t'); group_A_start.append(tmpMisure[0])
					HAM_startFILE_FID.write(str(tmpMisure[1]) + '\t'); group_HAM_start.append(tmpMisure[1])
					EUC_startFILE_FID.write(str(tmpMisure[2]) + '\t'); group_EUC_start.append(tmpMisure[2])	
					# PREVIOUS ONE (NO INFLUX) Defining concentration of the two vectors
					tmpMisure = distanceMisures(seqNOINFLUX, conc, seqOLDNOINFLUX, concOLD, idS)
					previousNOINFLUX_FILE_FID.write(str(tmpMisure[0]) + '\t'); group_A_prev_NI.append(tmpMisure[0])
					HAM_previousNOINFLUX_FILE_FID.write(str(tmpMisure[1]) + '\t'); group_HAM_prev_NI.append(tmpMisure[1])
					EUC_previousNOINFLUX_FILE_FID.write(str(tmpMisure[2]) + '\t'); group_EUC_prev_NI.append(tmpMisure[2])		
					# START (NO INFLUX) Defining concentration of the two vectors
					tmpMisure = distanceMisures(seqNOINFLUX, conc, seqSTART_NOINFLUX, concSTART, idS)
					startNOINFLUX_FILE_FID.write(str(tmpMisure[0]) + '\t'); group_A_start_NI.append(tmpMisure[0])
					HAM_startNOINFLUX_FILE_FID.write(str(tmpMisure[1]) + '\t'); group_HAM_start_NI.append(tmpMisure[1])
					EUC_startNOINFLUX_FILE_FID.write(str(tmpMisure[2]) + '\t'); group_EUC_start_NI.append(tmpMisure[2])
					
					# Compute and save floating average values
					if idS % 10 == 0:
						if idS != 0:
							previousFILE_FID_group.write(str(np.mean(group_A_prev)) + '\t')
							previousNOINFLUX_FILE_FID_group.write(str(np.mean(group_A_prev_NI)) + '\t')
							startFILE_FID_group.write(str(np.mean(group_A_start)) + '\t')
							startNOINFLUX_FILE_FID_group.write(str(np.mean(group_A_start_NI)) + '\t')
							HAM_previousFILE_FID_group.write(str(np.mean(group_HAM_prev)) + '\t')
							HAM_previousNOINFLUX_FILE_FID_group.write(str(np.mean(group_HAM_prev_NI)) + '\t')
							HAM_startFILE_FID_group.write(str(np.mean(group_HAM_start)) + '\t')
							HAM_startNOINFLUX_FILE_FID_group.write(str(np.mean(group_HAM_start_NI)) + '\t')
							EUC_previousFILE_FID_group.write(str(np.mean(group_EUC_prev)) + '\t')
							EUC_previousNOINFLUX_FILE_FID_group.write(str(np.mean(group_EUC_prev_NI)) + '\t')
							EUC_startFILE_FID_group.write(str(np.mean(group_EUC_start)) + '\t')
							EUC_startNOINFLUX_FILE_FID_group.write(str(np.mean(group_EUC_start_NI)) + '\t')	
							group_A_prev = []; group_HAM_prev = []; group_EUC_prev = [];
							group_A_start = []; group_HAM_start = []; group_EUC_start = [];
							group_A_prev_NI = []; group_HAM_prev_NI = []; group_EUC_prev_NI = [];
							group_A_start_NI = []; group_HAM_start_NI = []; group_EUC_start_NI = [];
							
							  
					# the new lists becomes the old one
					seqOLD = seq[:]
					seqOLDNOINFLUX = seqNOINFLUX[:]
					concOLD = conc[:]
					# Close file species		  					  
				  	fidSpecies.close()
				  
				  # Write newline char on result files
				  previousFILE_FID.write('\n')
				  previousNOINFLUX_FILE_FID.write('\n')
				  startFILE_FID.write('\n')
				  startNOINFLUX_FILE_FID.write('\n')
				  newSpecies_FID.write('\n')
				  livingSpecies_FID.write('\n')
				  totMass_FID.write('\n')
				  evaluatedFID.write('\n')
				  HAM_previousFILE_FID.write('\n')
				  HAM_previousNOINFLUX_FILE_FID.write('\n')
				  HAM_startFILE_FID.write('\n')
				  HAM_startNOINFLUX_FILE_FID.write('\n')
				  EUC_previousFILE_FID.write('\n')
				  EUC_previousNOINFLUX_FILE_FID.write('\n')
				  EUC_startFILE_FID.write('\n')
				  EUC_startNOINFLUX_FILE_FID.write('\n')
				  previousFILE_FID_group.write('\n')
				  previousNOINFLUX_FILE_FID_group.write('\n')
				  startFILE_FID_group.write('\n')
				  startNOINFLUX_FILE_FID_group.write('\n')
				  HAM_previousFILE_FID_group.write('\n')
				  HAM_previousNOINFLUX_FILE_FID_group.write('\n')
				  HAM_startFILE_FID_group.write('\n')
				  HAM_startNOINFLUX_FILE_FID_group.write('\n')
				  EUC_previousFILE_FID_group.write('\n')
				  EUC_previousNOINFLUX_FILE_FID_group.write('\n')
				  EUC_startFILE_FID_group.write('\n')
				  EUC_startNOINFLUX_FILE_FID_group.write('\n')
				  # Write zeroOneList on file
				  for zol in zeroList:
					strtoW = str(zol) + '\t'
					zeroOneSpeciesFID.write(strtoW)
				  zeroOneSpeciesFID.write('\n')
		else: 
			print " |- no result folder has been found"

# CLOSE FILES
previousFILE_FID.close()
previousNOINFLUX_FILE_FID.close()
startFILE_FID.close()
startNOINFLUX_FILE_FID.close()

HAM_previousFILE_FID.close()
HAM_previousNOINFLUX_FILE_FID.close()
HAM_startFILE_FID.close()
HAM_startNOINFLUX_FILE_FID.close()

EUC_previousFILE_FID.close()
EUC_previousNOINFLUX_FILE_FID.close()
EUC_startFILE_FID.close()
EUC_startNOINFLUX_FILE_FID.close()

previousFILE_FID_group.close()
previousNOINFLUX_FILE_FID_group.close()
startFILE_FID_group.close()
startNOINFLUX_FILE_FID_group.close()

HAM_previousFILE_FID_group.close()
HAM_previousNOINFLUX_FILE_FID_group.close()
HAM_startFILE_FID_group.close()
HAM_startNOINFLUX_FILE_FID_group.close()

EUC_previousFILE_FID_group.close()
EUC_previousNOINFLUX_FILE_FID_group.close()
EUC_startFILE_FID_group.close()
EUC_startNOINFLUX_FILE_FID_group.close()

newSpecies_FID.close()
livingSpecies_FID.close()
zeroOneSpeciesFID.close()
totMass_FID.close()
evaluatedFID.close()

print '|- FINISHED... SEE YOU NEXT TIME'

					
					

			
				

	
	
