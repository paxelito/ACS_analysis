#!/usr/bin/env python
# -*- coding: latin-1 -*-

import sys, os # Standard librar
import glob
import linecache


def zeroBeforeStrNum(tmpl, tmpL):
	strZero = ''
	nZeros = len(str(tmpL)) - len(str(tmpl))
	if nZeros > 0:
		for i in range(0,nZeros): strZero = strZero + '0'
	return strZero
#eof

try:
	StrFrom = sys.argv[1] # Here the path of the simulation output files
	StrTo = sys.argv[2] # Here the path of the new simulation files 
except:
	print "Usage:",sys.argv[0], "infile outfile"; sys.exit(1)

# Create absolute paths
StrFrom = os.path.abspath(StrFrom)	
StrTo = os.path.abspath(StrTo)	
origin = os.getcwd()
_LASTSPECIES_ = 29 
_REVRCTS_ = 1
_RATIOREV_ = 1000
_CLEAVAGE_ = 25.0
_CONDENSATION_ = 50.0
_COMPLEXFORM_ = 50.0

#ÊGo to the source folder
os.chdir(strFrom)

# Go into the result folder 
sourceResFolder = os.path.join(strFrom,"res")
os.chdir(sourceResFolder)


# Select last species, reactions and catalysis file
lastSpeciesFile = sort(glob.glob('species_*'))
lastReactionsFile = sort(glob.glob('reactions_*'))
lastCatalysisFile = sort(glob.glob('catalysis_*'))

# Move files into the new folder
fileDest = os.path.join(StrTo,"_acsinflux.csv")
os.system ("cp %s %s" % ("_acsinflux.csv",fileDest));
fileDest = os.path.join(StrTo,"_acsnrgbooleanfunctions.csv")
os.system ("cp %s %s" % ("_acsnrgbooleanfunctions.csv",fileDest));
fileDest = os.path.join(StrTo,"acsm2s.conf")
os.system ("cp %s %s" % ("acsm2s.conf",fileDest));

fileDest = os.path.join(StrTo,"_acsspecies.csv")
os.system ("cp %s %s" % (lastSpeciesFile,fileDest));

fileDest = os.path.join(StrTo,"_acsreactions.csv")
os.system ("cp %s %s" % (lastReactionsFile,fileDest));

fileDest = os.path.join(StrTo,"_acscatalysis.csv")
os.system ("cp %s %s" % (lastCatalysisFile,fileDest));

#ÊGo into the new folder
#ÊGo to the source folder
os.chdir(strTo)

# Reset Config File
# ACSCONF FILE
mod = open("acsm2s.conf").readlines()	
id = 0
for line in mod:
	if line[0] <> "#":
		linesplitted = line.split("=")
		if linesplitted[0] == 'randomSeed':
			linesplitted[1] = '0'
		if _REVRCTS_ == 1:
			if linesplitted[0] == 'reverseReactions':
				linesplitted[1] = '1'
			if linesplitted[0] == 'revRctRatio':
				linesplitted[1] = str(_RATIOREV_)		
		mod[id] = "=".join(linesplitted)
	id += 1	 
	
try:
	file = open("acsm2s.conf", "w")
	file.writelines(mod)
	file.close()
except IOError:
	print "Couldn't save configuration file"	

# Reset Species File 
mod = open("_acsspecies.csv").readlines()	
id = 0
for line in mod:
	linesplitted = line.split("\t")
	# Set concentration
	if id > _LASTSPECIES_:
		linesplitted[2] = '0'
	else:
		linesplitted[2] = '0.00110924'
	linesplitted[8] = '0' # age
	linesplitted[9] = '0' # reborn
	mod[id] = "\t".join(linesplitted)
	id += 1	 
try:
	file = open("_acsspecies.csv", "w")
	file.writelines(mod)
	file.close()
except IOError:
	print "Couldn't save species file"
	
#reset reactions file
mod = open("_acsreactions.csv").readlines()	
id = 0
for line in mod:
	linesplitted = line.split("\t")
	linesplitted[5] = '0' # counter
	mod[id] = "\t".join(linesplitted)
	id += 1	 
try:
	file = open("_acsreactions.csv", "w")
	file.writelines(mod)
	file.close()
except IOError:
	print "Couldn't save reactions file"
	
# reset and fix catalysis file
mod = open("_acscatalysis.csv").readlines()	
mod_rct = open("_acsreactions.csv").readlines()
id = 0
for line in mod_cat:
	flag = 0
	linesplitted = line.split("\t")
	linesplitted[3] = '0' # counter
	
	# Extract reaction from reaction file
	catRct = linecache.getline('_acsreactions.csv', int(linesplitted[2]))
	carRctSplit = catRct.split("\t")
	
	# Check whether condensation or cleavage
	if carRctSplit[1] == 0:
		linesplitted[4] = str(_CONDENSATION_) # K_cond
		linesplitted[5] = str(_CLEAVAGE_ / _RATIOREV_) # K_cleavage
		linesplitted[6] = str(_COMPLEXFORM_) + '\n' # K_complex
	else:
		linesplitted[4] = str(_CONDENSATION_ / _RATIOREV_) # K_cond
		linesplitted[5] = str(_CLEAVAGE_) # K_cleavage
		linesplitted[6] = str(_COMPLEXFORM_ / _RATIOREV_) + '\n' # K_complex

	mod[id] = "\t".join(linesplitted)
	id += 1	 
try:
	file = open("_acscatalysis.csv", "w")
	file.writelines(mod)
	file.close()
except IOError:
	print "Couldn't save catalysis file"


# Come back to the original folder
os.chdir(origin)

