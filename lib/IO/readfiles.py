#!/usr/bin/python
# -*- coding: latin-1 -*-

import sys, os # Standard librar
import glob
import numpy as np # Scientific library
from numpy import * 
from argparse import ArgumentParser
try:
    from pylab import *
except:
    pass
   
def readConfFile(tmpPath):
	#Open Parameter File
	paramFile = os.path.join(tmpPath,"acsm2s.conf")
	try:
		fid = open(paramFile, 'r')
	except:
		print 'impossible to load ', paramFile; sys.exit(1)
					
	# Read reaction probability from configuration file
	for line in fid:
		strLine = line.split('=')
		if strLine[0] == "nGEN":
			ngen = int(strLine[1])
		if strLine[0] == "nSIM":
			nsim = int(strLine[1])
		if strLine[0] == "nSeconds":
			totTimes = int(strLine[1])
		if strLine[0] == "energy":
			nrgType = int(strLine[1])
		if strLine[0] == "nReactions":
			totalRcts = int(strLine[1])
		if strLine[0] == "ECConcentration":
			nrgConc = float(strLine[1])
		if strLine[0] == "influx_rate":
			influx_rate = float(strLine[1])
		if strLine[0] == "maxLOut":
			maxLOut = int(strLine[1])
		if strLine[0] == "lastFiringDiskSpeciesID":
			lfdsID = int(strLine[1])
		if strLine[0] == "reactionProbability":
			rctProb = float(strLine[1])	
		if strLine[0] == "newSpeciesProbMinThreshold":
			newSpeciesProbMinThreshold = float(strLine[1])	
			
	return (ngen,nsim,totTimes,nrgType,totalRcts,nrgConc,influx_rate,maxLOut,lfdsID,rctProb,newSpeciesProbMinThreshold)

def read_sims_conf_file(paramFile):
	#Open Parameter File
	paramFile = os.path.abspath(paramFile)

	try:
		fid = open(paramFile, 'r')
	except:
		print 'impossible to load ', paramFile; sys.exit(1)
					
	# Read reaction probability from configuration file
	for line in fid:
		strLine = line.split('=')
		if strLine[0] == "nGEN":
			nGEN = int(strLine[1])
		if strLine[0] == "nSIM":
			nSIM = int(strLine[1])
		if strLine[0] == "nSeconds":
			nSeconds = int(strLine[1])
		if strLine[0] == "nReactions":
			nReactions = int(strLine[1])
		if strLine[0] == "nHours":
			nHours = int(strLine[1])
		if strLine[0] == "nAttempts":
			nAttempts = int(strLine[1])
		if strLine[0] == "randomSeed":
			randomSeed = int(strLine[1])
		if strLine[0] == "debugLevel":
			debugLevel = int(strLine[1])
		if strLine[0] == "timeStructuresSavingInterval":
			timeStructuresSavingInterval = float(strLine[1])
		if strLine[0] == "fileTimesSaveInterval":
			fileTimesSaveInterval = float(strLine[1])	
		if strLine[0] == "newSpeciesProbMinThreshold":
			newSpeciesProbMinThreshold = float(strLine[1])
		if strLine[0] == "lastFiringDiskSpeciesID":
			lastFiringDiskSpeciesID = int(strLine[1])	
		if strLine[0] == "overallConcentration":
			overallConcentration = float(strLine[1])
		if strLine[0] == "ECConcentration":
			ECConcentration = float(strLine[1])
		if strLine[0] == "alphabet":
			alphabet = str(strLine[1][:-1])
		if strLine[0] == "volume":
			volume = float(strLine[1])	
		if strLine[0] == "volumeGrowth":
			volumeGrowth = int(strLine[1])
		if strLine[0] == "stochDivision":
			stochDivision = int(strLine[1])
		if strLine[0] == "energy":
			energy = int(strLine[1])
		if strLine[0] == "ratioSpeciesEnergizable":
			ratioSpeciesEnergizable = float(strLine[1])
		if strLine[0] == "nonCatalyticMaxLength":
			nonCatalyticMaxLength = int(strLine[1])
		if strLine[0] == "reactionProbability":
			reactionProbability = float(strLine[1])	
		if strLine[0] == "cleavageProbability":
			cleavageProbability = float(strLine[1])
		if strLine[0] == "main_rev_rct_allowed":
			main_rev_rct_allowed = int(strLine[1])
		if strLine[0] == "reverseReactions":
			reverseReactions = int(strLine[1])	
		if strLine[0] == "revRctRatio":
			revRctRatio = float(strLine[1])
		if strLine[0] == "spontRct":
			spontRct = float(strLine[1])
		if strLine[0] == "K_ass":
			K_ass = float(strLine[1])
		if strLine[0] == "K_diss":
			K_diss = float(strLine[1])	
		if strLine[0] == "K_cpx":
			K_cpx = float(strLine[1])
		if strLine[0] == "K_cpxDiss":
			K_cpxDiss = float(strLine[1])
		if strLine[0] == "K_nrg":
			K_nrg = float(strLine[1])
		if strLine[0] == "K_nrg_decay":
			K_nrg_decay = float(strLine[1])
		if strLine[0] == "K_spont_ass":
			K_spont_ass = float(strLine[1])
		if strLine[0] == "K_spont_diss":
			K_spont_diss = float(strLine[1])	
		if strLine[0] == "moleculeDecay_KineticConstant":
			moleculeDecay_KineticConstant = float(strLine[1])
		if strLine[0] == "diffusion_contribute":
			diffusion_contribute = float(strLine[1])
		if strLine[0] == "solubility_threshold":
			solubility_threshold = float(strLine[1])		
		if strLine[0] == "influx_rate":
			influx_rate = float(strLine[1])
		if strLine[0] == "maxLOut":
			maxLOut = int(strLine[1])				
			
	#r1 = 0 : 9 || 10 : 17 || 11 : 24 || 25 : 36 || 37 : 39
	return (nGEN,nSIM,nSeconds,nReactions,nHours,nAttempts,randomSeed,debugLevel,timeStructuresSavingInterval,fileTimesSaveInterval,\
		    newSpeciesProbMinThreshold,lastFiringDiskSpeciesID,overallConcentration,ECConcentration,alphabet,volume,volumeGrowth,stochDivision,\
		    energy,ratioSpeciesEnergizable,nonCatalyticMaxLength,reactionProbability,cleavageProbability,main_rev_rct_allowed,reverseReactions,\
		    revRctRatio,spontRct,K_ass,K_diss,K_cpx,K_cpxDiss,K_nrg,K_nrg_decay,K_spont_ass,K_spont_diss,moleculeDecay_KineticConstant,\
		    diffusion_contribute,solubility_threshold,influx_rate,maxLOut)

def readInitConfFile(tmpPath):
	#Open Parameter File
	paramFile = os.path.join(tmpPath,"acsm2s.conf")
	try:
		fid = open(paramFile, 'r')
	except:
		print 'impossible to load ', paramFile; sys.exit(1)
					
	# Read reaction probability from configuration file
	for line in fid:
		strLine = line.split('=')
		if strLine[0] == "nGEN":
			ngen = int(strLine[1])
		if strLine[0] == "nSIM":
			nsim = int(strLine[1])
		if strLine[0] == "nSeconds":
			totTimes = int(strLine[1])
		if strLine[0] == "energy":
			nrgType = int(strLine[1])
		if strLine[0] == "nReactions":
			totalRcts = int(strLine[1])
		if strLine[0] == "ECConcentration":
			nrgConc = float(strLine[1])
		if strLine[0] == "influx_rate":
			influx_rate = float(strLine[1])
		if strLine[0] == "maxLOut":
			maxLOut = int(strLine[1])
		if strLine[0] == "lastFiringDiskSpeciesID":
			lfdsID = int(strLine[1])
		if strLine[0] == "reactionProbability":
			rctProb = float(strLine[1])	
			
	return (ngen,nsim,totTimes,nrgType,totalRcts,nrgConc,influx_rate,maxLOut,lfdsID,rctProb)

# Return Buffered species IDs
def readBufferedID(tmpPath):
	speciesFile = os.path.join(tmpPath,"_acsspecies.csv")
	try:
		fid = open(speciesFile, 'r')
	except:
		print 'impossible to load ', speciesFile; sys.exit(1)
		
	tempBufIDs = []
	for sp in fid:
		tmpID, tmpSeq, tmpConc, tmpDiff, tmpSol, tmpCpxDiss, tmpCpxCut, tmpEval, tmpAge, tmpReb, tmpCatID, tmpSubID, tmpKpho, tmpLoadConc, tmpConcLock = sp.split()
		if int(tmpConcLock) == 1:
			tempBufIDs.append(int(tmpID))
			
	return tempBufIDs

# Return CSTR flux
def readCSTRflux(tmpPath): 
	influxFile = os.path.join(tmpPath,"_acsinflux.csv")
	try:
		fid = open(influxFile, 'r')
	except:
		print 'impossible to load ', influxFile; sys.exit(1)
		
	tempIDs = []
	for sp in fid:
		tmpID, tmpProb = sp.split()
		tempIDs.append(int(tmpID))
			
	return tempIDs

# Return Reaction or Catalysis file
def loadAllData(tmpPath, tmpFname):
	
	fileName = os.path.join(tmpPath,tmpFname)
	data = np.loadtxt(fileName, dtype=float)
	return data

# Fixing the name of the file to create the generation number according to the total number of generations
def zeroBeforeStrNum(tmpl, tmpL):
	strZero = ''
	nZeros = len(str(tmpL)) - len(str(tmpl))
	if nZeros > 0:
		for i in range(0,nZeros): strZero = strZero + '0'
	return strZero

# function to split the single parameter reaction files in their compoenents
def splitRctParsLine(tmpLine):
	tmpReaction, tmpTime, tmpRtype, tmpCat, tmpMol_I, tmpMol_II, tmpMol_III, tmpLoadedMols,\
	tmpLoadedMolsConc, tmpGillMean, tmpGillSD, tmpGillEntropy, tmpNSCprob, tmpRevProb = tmpLine.split()
	
	
	if (int(tmpRtype) == 10) | (int(tmpRtype) == 11):
		tempVar = tmpMol_III
		tmpMol_III = tmpMol_I
		tmpMol_I = tmpCat
		tmpMol_II = tempVar
		
		
	
	return int(tmpReaction), float(tmpTime), int(tmpRtype), int(tmpCat), int(tmpMol_I), int(tmpMol_II), int(tmpMol_III), float(tmpLoadedMolsConc),\
		    int(tmpLoadedMols), float(tmpGillMean), float(tmpGillSD),  float(tmpGillEntropy), float(tmpNSCprob), float(tmpRevProb)


