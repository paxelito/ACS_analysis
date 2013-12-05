	#!/usr/bin/env python
# -*- coding: latin-1 -*-
'''Script to order the analysis of the divergences in time. 
'''
import sys, os # Standard librar
import glob
from copy import deepcopy
import numpy as np # Scientific library
import itertools as it
import random as ran

from time import time
from numpy import * 
from argparse import ArgumentParser
try:
    from pylab import *
except:
    pass
   
from lib.graph import raf
from lib.graph import scc
from lib.graph import network
from lib.dyn import dynamics as dn
from lib.model import species as sp
from lib.model import reactions
from lib.IO import *


#�Input parameters definition 
if __name__ == '__main__':
	parser = ArgumentParser(
				description='This script re-arrange results in a more friendly way from the angle analysis in time.'
				, epilog='''File with angle trajectories are created. ''') 
	parser.add_argument('-a', '--sysType', help='System Architecture [1:CLOSE, 2:PROTO, 3:CSTR], deafult: 1', default='1')
	parser.add_argument('-k', '--creationMethod', help='Network creation method (1: Filisetti, 2: Wim, 3: WimNoRevs, DEF: 1)', default='1', type=int)
	parser.add_argument('-f', '--lastFood', type=int, help='max food species length (deafult: 2)', default='2')
	parser.add_argument('-n', '--noCat', help='Non catalytic max length (default: 2)', default='5', type=int)
	parser.add_argument('-N', '--initAmount', help='Default Initial Amount', default='600', type=int)
	parser.add_argument('-s', '--initSet', type=int, help='Dimension of the initial set (Default: 4)', default='4')
	parser.add_argument('-m', '--maxDim', help='Max Dimension of the systems (Default: 6)', default='6', type=int)
	parser.add_argument('-o', '--strOut', help='Path for output file storing (Default: ./)', default='./')
	parser.add_argument('-I', '--conf', help='Configuration File (Default: ./acsm2s.conf)', default='./acsm2s.conf')
	parser.add_argument('-H', '--chemistries', help='Number of distinct chemistries to create', type=int, default='4')
	parser.add_argument('-i', '--iteration', help='Number of initial conditions (Default: 1)', default='1', type=int)
 	parser.add_argument('-v', '--avgCon', help='Catalysis level (deafult: 1), i.e. average reactions catalyzed per species', type=float, default='1')
	parser.add_argument('-c', '--rctRatio', help='Ratio between cleavages and condensations (default: 0, it means that the actual ratio is used)', default='0', type=float)
	parser.add_argument('-C', '--core', help='Number of core on which simulations are distributed', default='2', type=int)	
	parser.add_argument('-F', '--folderName', help='Simulation Folder Name (Deafault: SIMS)', default='SIMS')
	parser.add_argument('-R', '--revRcts', help='Reverse reactions are allowed to be created for chance(1: Yes, 0: No, Deafult: No)', default='0', type=int)	
	parser.add_argument('-r', '--randomSeed', help='random seed', type=int, default=None)
	args = parser.parse_args()
	
	_AVOGADRO_ = 6.022141e23
	
	# SIMULATION FOLDER CREATION
	folderName = os.path.join(args.strOut, args.folderName)
	if not os.path.isdir(folderName):
		try:
			os.mkdir(folderName)
		except:
			print "Impossible to create statistic directory", folderName; sys.exit(1)
			
	#�FOR EACH NETWORK
	fid_run = []
	for sngCore in range(args.core):
		zeroBeforeName =  readfiles.zeroBeforeStrNum(sngCore+1, args.core)
		runFileName = zeroBeforeName + str(sngCore+1) + '_simulation.sh'
		fname_run = os.path.join(folderName, runFileName)
		fid_run.append(open(fname_run, 'w'))
		
	# Read Conf File
	parameters = readfiles.read_sims_conf_file(args.conf)
	
	fidid = 0 # Core on which the chemistry will run
	for singlechem in range(args.chemistries):
		# CREATE CHEMIST FOLDER
		print "|- ARTIFICIAL CHEMISTRY ", singlechem+1, " Creation..."
		zeroBeforeName =  readfiles.zeroBeforeStrNum(singlechem+1, args.chemistries)
		chemFolder = 'CH' + str(zeroBeforeName) + str(singlechem+1)
		chemFolderPath = os.path.join(folderName, chemFolder)
		if not os.path.isdir(chemFolderPath):
			try:
				os.mkdir(chemFolderPath)
			except:
				print "Impossible to create statistic directory", chemFolderPath; sys.exit(1)	
		
		# CREATE ARTIFICIAL CHEMISTRY
		# Create foodlist
		foodList = range(int(2**(args.lastFood+1)-2))
		# Create species List
		speciesList = sp.createCompleteSpeciesPopulation(args.maxDim, parameters[14])
		# Crate chemist
		# compute population cardinality
		totSpecies = sp.getTotNumberOfSpeciesFromCompletePop(args.maxDim)
		
		# Compute overall conceivable number of reactions
		totCleavage = reactions.getNumOfCleavages(speciesList)
		
		if args.creationMethod == 1: totCond = reactions.getNumOfCondensations(totSpecies)
		else: totCond = 0
		totRcts = totCleavage + totCond
		
		# If the reaction probability is 0, it is set to the critical value of 1 reaction per species on average
		if parameters[21] == 0: prob = (1 / float(totRcts)) * args.avgCon
		else: prob = parameters[21]
		rctToCat = int(round(totSpecies * args.avgCon))

		# Create chemistry 
		chemFound = False
		if singlechem < (args.chemistries / 2.0):
			scanned = 0
			while chemFound == False:
				if scanned % 100 == 0: print "\t\t tried chemistries to find RAFs -> ", scanned
				rcts, cats = network.create_chemistry(args, speciesList, parameters, rctToCat, totCleavage, totCond)
				food = deepcopy(foodList)
				rafset = raf.rafsearch(rcts, cats, food) #�RAF search
				if len(rafset[2]) > 0:
					print rafset
					chemFound = True
				scanned += 1
		else:
			scanned = 0
			while chemFound == False:
				if scanned % 100 == 0: print "\t\t tried chemistries -> ", scanned
				rcts, cats = network.create_chemistry(args, speciesList, parameters, rctToCat, totCleavage, totCond)
				food = deepcopy(foodList)
				rafset = raf.rafsearch(rcts, cats, food) #�RAF search
				if len(rafset[2]) == 0:
					chemFound = True
				scanned += 1
		
		# CREATE DISTINCT INITIAL CONDITIONS		
		for singleCond in range(args.iteration):
			print " :- Iteration ", singleCond+1, " Creation..."
			# CREATE CONDITION FOLDER
			zeroBeforeName =  readfiles.zeroBeforeStrNum(singleCond+1, args.iteration)
			tmpCondFolder = zeroBeforeName + str(singleCond+1) + '_condition'
		
			condFolderPath = os.path.join(chemFolderPath, tmpCondFolder)
			if not os.path.isdir(condFolderPath):
				try:
					os.mkdir(condFolderPath)
				except:
					print "Impossible to create statistic directory", condFolderPath; sys.exit(1)
					
			resFolder = os.path.join(condFolderPath, 'res')
			if not os.path.isdir(resFolder):
				try:
					os.mkdir(resFolder)
				except:
					print "Impossible to create statistic directory", resFolder; sys.exit(1)
			
			#�-----------------------------
			# ARTIFICIAL CHEMISTRY CREATION
			# -----------------------------
			writefiles.write_acsms_file(condFolderPath,*parameters) # Save config file. 
			writefiles.write_and_create_std_nrgFile(condFolderPath) #�Save energy file.
			sp.createFileSpecies(condFolderPath, args, parameters, singleCond, speciesList)
			writefiles.write_and_createInfluxFile(condFolderPath, foodList)		
			writefiles.write_acsCatalysis_file(condFolderPath, cats)	
			writefiles.write_acsReactions_file(condFolderPath, rcts)	
			
			# UPDATE SIMULATION LUNCHER
			str2w = "echo 'Simulation " + chemFolderPath + "'\nnice ./carness ." + condFolderPath[4:] + "/ ." + condFolderPath + "/res/ ." + condFolderPath +\
					"/ > " + chemFolder + '_' + tmpCondFolder + ".log\n"
			
			fid_run[fidid].write(str2w)
			fidid += 1
			if fidid == args.core: fidid = 0 
	
	# Close fid runs		
	map(close(),fid_run)