#!/usr/bin/env python
# -*- coding: latin-1 -*-

import sys, os # Standard librar
import glob
import linecache
from argparse import ArgumentParser


def zeroBeforeStrNum(tmpl, tmpL):
	strZero = ''
	nZeros = len(str(tmpL)) - len(str(tmpl))
	if nZeros > 0:
		for i in range(0,nZeros): strZero = strZero + '0'
	return strZero
#eof

# Input parameters definition 
if __name__ == '__main__':
	parser = ArgumentParser(
				description='Script to create new init simulation files starting from the end of a previous simulation.'
				, epilog='''Species concentration are initilized according to the sequence uploaded, pay attention to the file selected. ''') 
	parser.add_argument('-f', '--StrFrom', help='Path of simulation files.', default='./')
	parser.add_argument('-t', '--StrTo', help='Path where new simulation init file will be stored', default='./')
	parser.add_argument('-s','--FileSpeciesToGetConc', help='Species file where concentrations to use are stored', default='')
	args = parser.parse_args()

# Create absolute paths
StrFrom = os.path.abspath(args.StrFrom)	
StrTo = os.path.abspath(args.StrTo)
StrFileSpeciesToGetConc = os.path.abspath(args.FileSpeciesToGetConc)
	
print '\n-> Origin Folder: ', StrFrom
print '-> Dest Folder: ', StrTo
print '-> StrFileSpeciesToGetConc: ', StrFileSpeciesToGetConc, '\n' 

origin = os.getcwd()
_LASTSPECIES_ = 29 
_REVRCTS_ = 1
_RATIOREV_ = 1000
_CLEAVAGE_ = 25.0
_CONDENSATION_ = 50.0
_COMPLEXFORM_ = 50.0
_INITSPECIESCONC_ = 0.00110924

# Go to the source folder
os.chdir(StrFrom)

# Move files into the new folder
fileDest = os.path.join(StrTo,"_acsinflux.csv")
os.system ("cp %s %s" % ("_acsinflux.csv",fileDest));
fileDest = os.path.join(StrTo,"_acsnrgbooleanfunctions.csv")
os.system ("cp %s %s" % ("_acsnrgbooleanfunctions.csv",fileDest));

# Go into the result folder 
sourceResFolder = os.path.join(StrFrom,"res")
os.chdir(sourceResFolder)

# Select last species, reactions and catalysis file
lastSpeciesFile = sorted(glob.glob('species_*'))
# the substring is reaction_1 in order to avoid the reaction_parameters file...
lastReactionsFile = sorted(glob.glob('reactions_0*'))
lastCatalysisFile = sorted(glob.glob('catalysis_*'))

fileDest = os.path.join(StrTo,"acsm2s.conf")
os.system ("cp %s %s" % ("acsm2s.conf",fileDest));

print "-> files _acsinflux.csv, _acsnrgbooleanfunctions.csv and acsm2s.conf have been copied...\n"

fileDest = os.path.join(StrTo,"_acsspecies.csv")
os.system ("cp %s %s" % (lastSpeciesFile[-1],fileDest));

fileDest = os.path.join(StrTo,"_acsreactions.csv")
os.system ("cp %s %s" % (lastReactionsFile[-1],fileDest));

fileDest = os.path.join(StrTo,"_acscatalysis.csv")
os.system ("cp %s %s" % (lastCatalysisFile[-1],fileDest));

print "-> files\n\t", lastSpeciesFile[-1], "\n\t", lastReactionsFile[-1], "\n\t", lastCatalysisFile[-1], "\n"
print "\thave been copied and renamed into the new folder" 

#ÊGo into the new folder
#ÊGo to the source folder
os.chdir(StrTo)

# Reset Config File
# ACSCONF FILE
mod = open("acsm2s.conf").readlines()	
id = 0
for line in mod:
	if line[0] <> "#":
		linesplitted = line.split("=")
		if linesplitted[0] == 'randomSeed':
			linesplitted[1] = '0\n'
		if linesplitted[0] == 'K_diss':
			linesplitted[1] = str(_CLEAVAGE_) + '\n'
		if linesplitted[0] == 'K_ass':
			linesplitted[1] = str(_CONDENSATION_) + '\n'
		if linesplitted[0] == 'K_cpx':
			linesplitted[1] = str(_COMPLEXFORM_) + '\n'		
		if _REVRCTS_ == 1:
			if linesplitted[0] == 'reverseReactions':
				linesplitted[1] = '1\n'
			if linesplitted[0] == 'revRctRatio':
				linesplitted[1] = str(_RATIOREV_) + '\n'		
		mod[id] = "=".join(linesplitted)
	id += 1	 
	
try:
	file = open("acsm2s.conf", "w")
	file.writelines(mod)
	file.close()
except IOError:
	print "Couldn't save configuration file"	
	
print "-> file acsm2s.conf has been processed for the new simulation..."

# In concentrations come from specific file they are stored into a vector
if args.FileSpeciesToGetConc != '':
	print "\nLoading concentration from file ", StrFileSpeciesToGetConc
	concs = []
	specFileLines = open(StrFileSpeciesToGetConc).readlines()
	for specFileLine in specFileLines:
		linesplitted = specFileLine.split("\t")
		concs.append(linesplitted[2]) 
# Reset Species File 
mod = open("_acsspecies.csv").readlines()	
id = 0
for line in mod:
	linesplitted = line.split("\t")
	# Set concentration
	if args.FileSpeciesToGetConc == '':
		if id > _LASTSPECIES_:
			linesplitted[2] = '0'
		else:
			linesplitted[2] = str(_INITSPECIESCONC_)
	else:
		if len(concs) > id:
			linesplitted[2] = str(concs[id])
		else:
			linesplitted[2] = '0'
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
	
print "-> file _acsspecies.csv has been processed for the new simulation..."
	
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

print "-> file _acsreactions.csv has been processed for the new simulation..."
	
# reset and fix catalysis file
mod = open("_acscatalysis.csv").readlines()	
mod_rct = open("_acsreactions.csv").readlines()
id = 0
for line in mod:
	flag = 0
	linesplitted = line.split("\t")
	linesplitted[3] = '0' # counter
	
	# Extract reaction from reaction file
	catRct = linecache.getline('_acsreactions.csv', int(linesplitted[2])+1)
	carRctSplit = catRct.split("\t")
	
	# Check whether condensation or cleavage
	if int(carRctSplit[1]) == 0:
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
	
print "-> file _acscatalysis.csv has been processed for the new simulation..."

print "NEW SIMULATION READY TO START\n"


# Come back to the original folder
os.chdir(origin)

