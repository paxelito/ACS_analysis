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
		if strLine[0] == "reactionProbability":
			rp = float(strLine[1])
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
			
	return (ngen,nsim,totTimes,nrgType,totalRcts,nrgConc,influx_rate,maxLOut,lfdsID)

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
	print fileName
	data = np.loadtxt(fileName, dtype=float)
	return data

	
		

