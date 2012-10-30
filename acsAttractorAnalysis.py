#!/usr/bin/env python
'''Function to compute the different attractors emerging from a simulazione. The algorithm compares all the final states of the simulation computing 
	The differences between those
'''

import sys, os # Standard library
import datetime as dt
import linecache as lc
import glob

try:
    from pylab import *
except:
    pass

def orderDictionary(tmpDict):
	newList = {}
	for key in sorted(tmpDict.iterkeys()):
		newList[key] = tmpDict[key]
	return newList
	
def orderDictionary2(adict): 
	keys = adict.keys() 
	keys.sort() 
	return [adict[key] for key in keys]
	
# get path for placing simulation
try:
	StrPath = sys.argv[1] # Here the path of the simulation output file
except:
	print "Usage:",sys.argv[0], "infile outfile"; sys.exit(1)
	
print "Simulation Results Path: ", StrPath

today = dt.date.today()
	
tmpDirs = os.listdir(StrPath)

print tmpDirs
	
for tmpDir in tmpDirs:
	os.chdir(StrPath)
	totDirName = StrPath + '/' + tmpDir
	if os.path.isdir(totDirName):
		# Move to the directory 
		os.chdir(totDirName)
		resDirPath = os.path.join(totDirName, "res")
		if os.path.isdir(resDirPath):
			os.chdir(resDirPath)
		
			# Searching for files
			speciesFiles = sorted(glob.glob(os.path.join(resDirPath,'species_*')))
			speciesFile = speciesFiles[-1]
		
			print " "
			print "|- Species File: ", speciesFile
			# Open Catalysis File
			try:
				fidSpecies = open(speciesFile, 'r')
			except:
				print 'impossible to load ', speciesFile; sys.exit(1)
			
			# For each last species file
			seq = []; conc = []; speciesConc = []
			for sp in fidSpecies:
				tmpID, tmpSeq, tmpConc, tmpDiff, tmpSol, tmpCpxDiss, tmpCpxCut, tmpEval, tmpAge, tmpReb, tmpCatID, tmpSubID, tmpKpho, tmpLoadConc, tmpConcLock = sp.split()
				if (int(tmpCpxCut) == 0) & (float(tmpConc) > 0):
					seq.append(str(tmpSeq))
					conc.append(float(tmpConc))
					speciesConc.append([str(tmpSeq),float(tmpConc)])
					
			#print speciesConc			
			print "-------------"
			speciesConc.sort()
			print " "
			print speciesConc
			print " "

			
    			


					
					
			fidSpecies.close()	
			
				

	
	
