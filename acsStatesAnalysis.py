#!/usr/bin/env python
# -*- coding: latin-1 -*-
'''Script to compute the distance between different state of the same simulation. 
	comparison between t0 and t, t-1 and t adopting three different distance misure:
	angle, euclidian distance and hamming distance. 
	Moreover the script make an analysis of all the aggregative variables. 
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
		tmpCos = float(np.dot(vecX,vecY) / (np.linalg.norm(vecX) * np.linalg.norm(vecY)))	
		if tmpCos >= 1.0: tmpCos = 1.0		
		strtoW[0] = np.arccos(tmpCos)  * 360.00 / 6.28
		
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


if __name__ == '__main__':
	parser = ArgumentParser(
				description='Function to evaluate the activity of each species during the simulation, \
				catalyst substrate product or nothing. Moreover the script recognize all those molecules functioning as hub'
				, epilog='''File with angle trajectories are created. ''') 
	parser.add_argument('-p', '--strPath', help='Path where files are stored', default='./')
	parser.add_argument('-l', '--tmpMaxFluxL', type=int, help='Max length of the incoming flux', default=2)
	args = parser.parse_args()
	
	print "Simulation Results Path: ", args.StrPath
	
	today = dt.date.today()
	
	StrPath = os.path.abspath(args.StrPath)
	tmpMaxFluxL = args.tmpMaxFluxL
		
	tmpDirs = sort(os.listdir(StrPath))
	
	os.chdir(StrPath)
	 
	print ''
	print 'ACS STATES ANALYSER'
	print ''
	print '|- STEP 1. Creating common sorted species list...'
	
	#currentDir = StrPath.split("/")[-1]
	currentDir = ''
	ndn = currentDir + '_0_new_allStatResults'
	newdirAllResults = os.path.join(os.curdir, ndn)
	if not os.path.isdir(newdirAllResults):
		try:
			os.mkdir(newdirAllResults)
		except:
			print "Impossible to create statistic directory", newdirAllResults; sys.exit(1)
			
	os.chdir(newdirAllResults)
	
	# Open Files containing results 
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
	
	ANG_middlePreviousFILE_FID = open('STAT_ANG_t_middle_NOINFLUX.csv', 'w')
	HAM_middlePreviousFILE_FID = open('STAT_HAM_t_middle_NOINFLUX.csv', 'w')
	EUC_middlePreviousFILE_FID = open('STAT_EUC_t_middle_NOINFLUX.csv', 'w')
	
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
	
	newSpecies_FID = open('STAT_GENERAL_newSpecies.csv', 'w')
	livingSpecies_FID = open('STAT_GENERAL_livingSpecies.csv', 'w')
	mols_FID = open('STAT_GENERAL_mols.csv', 'w')
	totMass_FID = open('STAT_GENERAL_overallMass.csv', 'w')
	totOverallMass_FID = open('STAT_GENERAL_overallTotMass.csv', 'w')
	complex_FID = open('STAT_GENERAL_complex.csv', 'w')
	complexMols_FID = open('STAT_GENERAL_complexMols.csv', 'w')
	evaluatedFID = open('STAT_GENERAL_evaluated.csv', 'w')
	zeroOneSpeciesFID = open('STAT_GENERAL_zeroOneSpecies.csv', 'w')
	biodeversityFID = open('STAT_GENERAL_bioDiversity.csv', 'w')
	
	os.chdir(StrPath)
	validDir = 1
	for IDdir, tmpDir in enumerate(tmpDirs):
	
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
				conf = readfiles.read_sims_conf_file() #ÊConfiguration file upload
				tmpVolume = conf[15]
				
				group_A_prev = []; group_HAM_prev = []; group_EUC_prev = [];
				group_A_start = []; group_HAM_start = []; group_EUC_start = [];
				group_A_prev_NI = []; group_HAM_prev_NI = []; group_EUC_prev_NI = [];
				group_A_start_NI = []; group_HAM_start_NI = []; group_EUC_start_NI = [];
				
				print "This Simulation is Composed of ", numberOfGen, " successive runs..."
				
				for ngen in range(1,numberOfGen+1):
				  
					  strZeros = zeroBeforeStrNum(ngen, numberOfGen)
					  
					  if ngen == 1:
					  	strSpeciesZero = 'species_' + strZeros + str(0) + '*'
					  	speciesFilesZero = sorted(glob.glob(os.path.join(resDirPath,strSpeciesZero)))
					  				  
					  strSpecies = 'species_' + strZeros + str(ngen) + '*'  
						  
					  # Searching for files
					  speciesFiles = sorted(glob.glob(os.path.join(resDirPath,strSpecies)))
					  
					  # Create a 0 list containing the present species
					  zeroList = returnZeroSpeciesList(speciesFiles[-1])
					  
					  if ngen == 1:
					  	speciesFiles = speciesFilesZero + speciesFiles
					  	
					  # Create matrix with dimension (NumberOfFileSpecies X Number of species (not complexes) )
					  # print len(speciesFiles)
					  speciesConcs = np.zeros((len(speciesFiles)+1,len(zeroList)),dtype=np.float64)
					  	
					  # Initialize moving average lists
					  seqOLD = []; seqOLDNOINFLUX = []; concOLD = []
					  seqSTART = []; seqSTART_NOINFLUX = []; concSTART = []
					  totMass = []; obsSpecies = [];
					  seqMIDDLE_NOINFLUX = []; concMIDDLE = [];
					  
					  oldNumberOfSpecies = 0
					  
					  # FOR EACH FILE SPECIES ---------------------------
					  print len(speciesFiles), " Species files will be evaluated..."
					  for idS, sngSpeciesFile in enumerate(speciesFiles):
					  	
					  	#print '  |- Species File: ', sngSpeciesFile	
						try:
							fidSpecies = open(sngSpeciesFile, 'r')
						except:
							print ' |- impossible to load ', sngSpeciesFile; sys.exit(1)
							
						seq = []; conc = []; seqNOINFLUX = []; numberOfSpecies = 0; tmpMass = 0; tmpObsSpecies = 0
						tmpMols = 0; tmpCpx = 0; tmpCpxMols = 0; tmpTotMass = 0
						for sp in fidSpecies:
							tmpID, tmpSeq, tmpConc, tmpDiff, tmpSol, tmpCpxDiss, tmpCpxCut, tmpEval, tmpAge, tmpReb, tmpCatID, tmpSubID, tmpKpho, tmpLoadConc, tmpConcLock = sp.split()
							if (int(tmpCpxCut) == 0) & (np.float64(tmpConc) > 0):
								seq.append(str(tmpSeq))
								conc.append(np.float64(tmpConc))
								if idS == 0:
									seqSTART.append(str(tmpSeq))
									concSTART.append(np.float64(tmpConc))
								if idS == 50:
									concMIDDLE.append(np.float64(tmpConc))							
								if len(str(tmpSeq)) > tmpMaxFluxL:
									seqNOINFLUX.append(str(tmpSeq))	
									if idS == 0:
										seqSTART_NOINFLUX.append(str(tmpSeq))
									if idS == 50:
										seqMIDDLE_NOINFLUX.append(str(tmpSeq))
								# Set 1 in zeroList (species has been created at least once)
								zeroList[numberOfSpecies] = 1
								# Update species concentrations matrix
								if idS == 0: speciesConcs[idS, numberOfSpecies] = tmpID
								if speciesConcs[0, numberOfSpecies] == 0: speciesConcs[0, numberOfSpecies] = tmpID
								speciesConcs[idS+1, numberOfSpecies] = tmpConc
								# Update systems mass
								tmpMass += len(str(tmpSeq)) * int(round(float(tmpConc) * 6.022e23 * tmpVolume))
								tmpMols += int(round(float(tmpConc) * 6.022e23 * tmpVolume))
							if (int(tmpCpxCut) > 0) & (float(tmpConc) > 0):
								tmpCpx += 1
								tmpCpxMols += int(round(float(tmpConc) * 6.022e23 * tmpVolume))
							if (float(tmpConc) > 0):
								tmpTotMass += len(str(tmpSeq)) * int(round(float(tmpConc) * 6.022e23 * tmpVolume))
							# If the species is not a complex, the number of species is updated	
							if int(tmpCpxCut) == 0:
								numberOfSpecies += 1	
							if int(tmpCpxCut) == 0 and int(tmpEval) == 1:
								tmpObsSpecies += 1		
						
						# Compute biodiversity index
						bioDivInd = 0
						for sngConc in conc:
							bioDivInd += pow(float(round(float(sngConc) * 6.022e23 * tmpVolume) / tmpMols),2)
						bioDivInd = float(1 / bioDivInd)
						
						# Compute number of new species		
						deltaNspecies = numberOfSpecies - oldNumberOfSpecies
						oldNumberOfSpecies = numberOfSpecies 
						strtoW = str(deltaNspecies) + '\t'
						newSpecies_FID.write(strtoW)	
						strtoW = str(len(conc)) + '\t'
						livingSpecies_FID.write(strtoW)		
						strtoW = str(tmpMass) + '\t'
						totMass_FID.write(strtoW)	
						strtoW = str(tmpMols) + '\t'
						mols_FID.write(strtoW)
						strtoW = str(tmpObsSpecies) + '\t'
						evaluatedFID.write(strtoW)	
						strtoW = str(bioDivInd) + '\t'
						biodeversityFID.write(strtoW)
						strtoW = str(tmpCpx) + '\t'
						complex_FID.write(strtoW)
						strtoW = str(tmpCpxMols) + '\t'
						complexMols_FID.write(strtoW)
						strtoW = str(tmpTotMass) + '\t'
						totOverallMass_FID.write(strtoW)			
						
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
						# MIDDLE (NO INFLUX) Defining concentration of the two vectors
						if idS >= 50:
							tmpMisure = distanceMisures(seqNOINFLUX, conc, seqMIDDLE_NOINFLUX, concMIDDLE, idS)
							ANG_middlePreviousFILE_FID.write(str(tmpMisure[0]) + '\t')
							HAM_middlePreviousFILE_FID.write(str(tmpMisure[1]) + '\t')
					  		EUC_middlePreviousFILE_FID.write(str(tmpMisure[2]) + '\t')
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
					  mols_FID.write('\n')
					  evaluatedFID.write('\n')
					  biodeversityFID.write('\n')
					  complex_FID.write('\n')
					  complexMols_FID.write('\n')
					  totOverallMass_FID.write('\n')
					  
					  HAM_previousFILE_FID.write('\n')
					  HAM_previousNOINFLUX_FILE_FID.write('\n')
					  HAM_startFILE_FID.write('\n')
					  HAM_startNOINFLUX_FILE_FID.write('\n')
					  
					  EUC_previousFILE_FID.write('\n')
					  EUC_previousNOINFLUX_FILE_FID.write('\n')
					  EUC_startFILE_FID.write('\n')
					  EUC_startNOINFLUX_FILE_FID.write('\n')
					  
					  if idS >= 50:
					  	ANG_middlePreviousFILE_FID.write('\n')
					  	HAM_middlePreviousFILE_FID.write('\n')
					  	EUC_middlePreviousFILE_FID.write('\n')
					  
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
					  
					  # clean concentrations matrix from always zero species
					  speciesConcs = speciesConcs[:,speciesConcs.sum(0)!=0]
					  # Save concentrations matrix
					  filename = "STAT_species_Concentrations_"  + str(validDir) + "_" + str(ngen) + ".csv"
					  np.savetxt(os.path.join(StrPath,ndn,filename), speciesConcs, delimiter='\t', fmt='%.10f')
				validDir += 1
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
	
	ANG_middlePreviousFILE_FID.close()
	HAM_middlePreviousFILE_FID.close()
	EUC_middlePreviousFILE_FID.close()
	
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
	mols_FID.close()
	biodeversityFID.close()
	complex_FID.close()
	complexMols_FID.close()
	totOverallMass_FID.close()
	
	print '|- FINISHED... SEE YOU NEXT TIME'

					
					

			
				

	
	
