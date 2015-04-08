	#!/usr/bin/env python
# -*- coding: latin-1 -*-
'''Script to initialize random catalytic nets 
python <path>/GIT/ACS_analysis/initializator.py -t2 -a0 -o ~/Documents/lavoro/protocell/init/ 
												-k3 -d2 -K10 -f2 -n2 -s6 -m6 -p7 -I ~/Documents/lavoro/protocell/init/acsm2s.conf 
												-N600 -B600 -x1 -O0 -H10 -v2.5 -c0.5 -F TEST -i 1 -S2 -u -P2 -S2 -A0.1
2015 - experiments on synchronization, init command
python <path>/initializator.py t2 -H20 -K20 -u -v1.0 -P2 -S2 -F <folder_name> -A0.01
and
python <path>/initializator.py t2 -H20 -K20 -u -v2.5 -P2 -S2 -F <folder_name> -A0.01
'''
import sys, os # Standard librar
#import glob
from copy import deepcopy
import numpy as np # Scientific library
#import itertools as it
import random as ran

#from time import time
#from numpy import * 
from argparse import ArgumentParser
#try:
#    from pylab import *
#except:
#    pass

from lib.graph import raf
from lib.graph import scc
from lib.graph import network
#from lib.dyn import dynamics as dn
from lib.model import species as sp
from lib.model import reactions
from lib.IO import *


#æInput parameters definition 
if __name__ == '__main__':
	parser = ArgumentParser(
				description='This script initialize new simulation structures.'
				, epilog='''File with angle trajectories are created. ''') 
	parser.add_argument('-F', '--folderName', help='Simulation Folder Name (Deafault: SIMS)', default='SIMS')
	parser.add_argument('-C', '--core', help='Number of core on which simulations are distributed (def:2)', default='2', type=int)
	parser.add_argument('-t', '--sysType', help='System Architecture [1:CLOSE, 2:PROTO, 3:CSTR], deafult: 1', default='2')
	parser.add_argument('-a', '--prefAttach', help='Type of catalyst choice (1: Preferential Attachment, 0: Random attachment, DEF: 0', default='0', type=int)
	parser.add_argument('-o', '--strOut', help='Path for output file storing (Default: ./)', default='./')
	parser.add_argument('-k', '--creationMethod', help='Network creation method (1: Filisetti, 2: Wim, 3: WimNoRevs, 4: WIM_RAFinREV_noRAFinNOrev, DEF: 3)', default='3', type=int)
	parser.add_argument('-R', '--revRcts', help='Reverse reactions are allowed to be created for chance(1: Yes, 0: No, Deafult: No)', default='0', type=int)		
	parser.add_argument('-d', '--directRctDirection', help='Direction of the forward reaction where necessary (1: cleavage, 0: condensation, 2: random with probability 0.5,Default: 2)', default='2', type=int)
	parser.add_argument('-H', '--chemistries', help='Number of distinct chemistries to create (Def: 4)', type=int, default='4')	
	parser.add_argument('-K', '--chemistriesWithRAF', help='Number of Chemistries with RAF to create (def:0), -1 does not care about RAFs', default='0', type=int)
	parser.add_argument('-u', '--autocat', help='Allow autocatalysis in principle (no param:TRUE, -u=FALSE)', action="store_false", default=True)
	parser.add_argument('-f', '--lastFood', type=int, help='max food species length (deafult: 2)', default='2')
	parser.add_argument('-s', '--initSet', type=int, help='Max Dimension of the initial set (Default: 6)', default='6')
	parser.add_argument('-m', '--maxDim', help='Max Dimension of the systems (Default: 6)', default='6', type=int)
	parser.add_argument('-n', '--noCat', help='Non catalytic max length (default: 2)', default='2', type=int)
	parser.add_argument('-p', '--redConc', help='Minimal dimension with reduced concentration. If it is greater than --initset no reduced concentration will be adopted (default: 7)', default='7', type=int)
	parser.add_argument('-N', '--initAmount', help='Default Initial Amount (def:600)', default='600', type=int)
	parser.add_argument('-B', '--initBufferAmount', help='Default Initial Buffer Amount (def:600)', default='600', type=int)
	parser.add_argument('-x', '--fixedConcentration', help='--initAmount is the average amount (0) or the real amount (1)  (def:0)', default='0', type=int)
	parser.add_argument('-O', '--holesperc', help='Percentage of holes in initial concentrations (def:0)', default='0', type=float)
	parser.add_argument('-I', '--conf', help='Configuration File (Default: ./acsm2s.conf)', default='./acsm2s.conf')
	parser.add_argument('-i', '--iteration', help='Number of initial conditions (Default: 1)', default='1', type=int)
	parser.add_argument('-v', '--avgCon', help='Catalysis level (deafult: 1), i.e. average reactions catalyzed per species', type=float, default='1')
	parser.add_argument('-c', '--rctRatio', help='Ratio between cleavages and condensations (default: 0.5)', default='0.5', type=float)
	parser.add_argument('-A', '--alpha', help='Kinetic rate of membrane growth by means of the catalytic activity of the catalytic molecoles (def:1)', default='1', type=float)
	parser.add_argument('-P', '--rafPresence', help='Force the presence of RAF of 1 dimension, i.e. self-catalysis (1), or bigger (2...N) or no RAF at all (0), Default: 0', default='0', type=int)
	parser.add_argument('-S', '--sccinraf', help='minimal dimension of the SCC within a RAF (def: 0)', type=int, default=0)
	parser.add_argument('-r', '--randomSeed', help='random seed', type=int, default=None)
	args = parser.parse_args()
	
	_AVOGADRO_ = 6.02214179e+23
	
	# SIMULATION FOLDER CREATION
	folderName = os.path.join(args.strOut, args.folderName)
	if not os.path.isdir(folderName):
		try:
			os.mkdir(folderName)
		except:
			print "Impossible to create statistic directory", folderName; sys.exit(1)
			
	# FOR EACH NETWORK
	fid_run = []
	for sngCore in range(args.core):
		zeroBeforeName =  readfiles.zeroBeforeStrNum(sngCore+1, args.core)
		runFileName = zeroBeforeName + str(sngCore+1) + '_simulation.sh'
		fname_run = os.path.join(folderName, runFileName)
		fid_run.append(open(fname_run, 'w'))
		
	# Read Conf File !!!!
	parameters = readfiles.read_sims_conf_file(args.conf)
	
	fidid = 0 # Core on which the chemistry will run
	
	fname_initRafRes = os.path.join(folderName, '0_theoreticalRAFanalysis.csv')
	fid_initRafRes = open(fname_initRafRes, 'w')
	strToWrite = "RAFsize\tSCCsize\tClosure\tCats\tRAF(unique)\tRAF\n"
	fid_initRafRes.write(strToWrite)
	
	for idchem, singlechem in enumerate(range(args.chemistries)):
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
		originalSpeciesList = deepcopy(speciesList)
		# Crate chemist
		# compute population cardinality
		totSpecies = sp.getTotNumberOfSpeciesFromCompletePop(args.maxDim)
		
		# Compute overall conceivable number of reactions
		totCleavage = reactions.getNumOfCleavages(originalSpeciesList)
		
		# Compute the number of reactions according to the chemistry type adopted
		if args.creationMethod == 1: totCond = reactions.getNumOfCondensations(totSpecies)
		else: totCond = 0
		totRcts = totCleavage + totCond
		
		# If the reaction probability is 0, it is set to the critical value of 1 reaction per species on average
		if parameters[21] == 0: prob = (1 / float(totRcts)) * args.avgCon
		else: prob = parameters[21]
		rctToCat = int(round(totSpecies * args.avgCon))

		# Create chemistry 
		chemFound = False
		if singlechem < args.chemistriesWithRAF:
			scanned = 0
			while chemFound == False:
				# CHEMISTRY CREATION
				if scanned % 100 == 0: print "\t\t tried chemistries to find RAFs -> ", scanned
				rcts, cats, speciesList, rcts_no_rev, cats_no_rev = network.create_chemistry(args, originalSpeciesList, parameters, rctToCat, totCleavage, totCond, rctToCat, args.autocat)
				food = deepcopy(foodList)
				rafset = raf.rafsearch(rcts, cats, food) # RAF search
				if args.creationMethod == 4: 
					food = deepcopy(foodList)
					rafset_no_rev = raf.rafsearch(rcts_no_rev, cats_no_rev, food) # RAF search
					
				# IF the RAF has been found the chemistry is valid, is creationMethod == 4 no RAF at all must be present in the chemistry without reverse reactions
				#if (len(rafset[2]) == args.rafPresence) & (len(set([6,8,13,16])&set(rafset[3]))==0):
				if (len(rafset[2]) >= args.rafPresence):
					if args.creationMethod == 4:
						if len(rafset_no_rev[2]) == 0:
							chemFound = True
					else:		
						chemFound = True
					# RAF SCC control 
					if chemFound == True:
						if len(rafset[2]) > 0: 
							rctsRAF = rcts[np.any(rcts[:, 0] == np.expand_dims(rafset[2],1), 0), :]
							catprodgraph = scc.createNetXGraphForRAF(rctsRAF, rafset[0], cats)
							if args.sccinraf >= 0:
								scc_in_raf = scc.checkMinimalSCCdimension(catprodgraph, args.sccinraf)
								#print scc_in_raf
								if scc_in_raf[0] == True:
									if args.sccinraf > 0:
										if args.sccinraf <= int(scc_in_raf[1][0]):
											scc.printSCConFile(scc_in_raf[2], folderName, idchem+1)
										else:
											chemFound = False 
											scc_in_raf = [0,0]
											raw_input("165")
									else:
										if sum(scc_in_raf[1]) == 0:
											scc.printSCConFile(scc_in_raf[2], folderName, idchem+1)
										else:
											chemFound = False 
											scc_in_raf = [0,0]
											raw_input("172")
								else:
									if args.sccinraf > 0:
										scc_in_raf = [0,0]
										chemFound = False
										raw_input("177")
									else:  
										scc_in_raf = [0,0]
							else:
								scc_in_raf = scc.checkMinimalSCCdimension(catprodgraph, args.sccinraf)
								#scc_in_raf = [0,0]
						else:
							scc_in_raf = [0,0]
						
				scanned += 1
		else:
			scanned = 0
			while chemFound == False:
				if scanned % 100 == 0: print "\t\t tried chemistries -> ", scanned
				rcts, cats, speciesList, rcts_norev, cats_norev = network.create_chemistry(args, originalSpeciesList, parameters, rctToCat, totCleavage, totCond, rctToCat, args.autocat)
				food = deepcopy(foodList)
				rafset = raf.rafsearch(rcts, cats, food) # RAF search
				if len(rafset[2]) > 0: # If the RAF contains at least one species
					rctsRAF = rcts[np.any(rcts[:, 0] == np.expand_dims(rafset[2],1), 0), :]
					catprodgraph = scc.createNetXGraphForRAF(rctsRAF, rafset[0], cats)
					if args.sccinraf >= 0:
						scc_in_raf = scc.checkMinimalSCCdimension(catprodgraph, args.sccinraf)
					
				if (len(rafset[2]) == 0) | (args.chemistriesWithRAF==-1):
					scc_in_raf = [0,0]
					chemFound = True
				scanned += 1
				
		strToWrite = str(len(rafset[2])) + "\t" + str(scc_in_raf[1]) + "\t" + str(rafset[0]) + "\t" + str(rafset[3]) + "\t" + str(rafset[4]) + "\t" + str(rafset[2]) + "\n"
		fid_initRafRes.write(strToWrite)
		
		# CREATE DISTINCT INITIAL CONDITIONS	
		# If volume growth define species contributing to the volume growth
		selcats = None
		if parameters[16] > 0: 
			if len(rafset[3]) > 0: 
				#selcats = [ran.choice(rafset[3]) for i in range(0,6)]
				selcats = rafset[3]
			else:
				selcats = ran.sample(range(2**(args.lastFood+1)-2,2**(args.maxDim+1)-2),6)
				
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
			
			# -----------------------------
			# ARTIFICIAL CHEMISTRY CREATION
			# -----------------------------
			writefiles.write_acsms_file(condFolderPath,*parameters) # Save config file. 
			writefiles.write_and_create_std_nrgFile(condFolderPath) # Save energy file.
			sp.createFileSpecies(condFolderPath, args, parameters, singleCond, speciesList, selcats)
			writefiles.write_and_createInfluxFile(condFolderPath, foodList)		
			writefiles.write_acsCatalysis_file(condFolderPath, cats)	
			writefiles.write_acsReactions_file(condFolderPath, rcts)	
			
			# UPDATE SIMULATION LUNCHER
			str2w = "echo 'Simulation " + condFolderPath + "'\nnice ./carness ." + condFolderPath + "/ ." + condFolderPath + "/res/ ." + condFolderPath +\
					"/ > " + chemFolder + '_' + tmpCondFolder + ".log\n"
			
			fid_run[fidid].write(str2w)
			fidid += 1
			if fidid == args.core: fidid = 0 
	
	# Close fid runs	
	fid_initRafRes.close()	
	map(lambda x: x.close(),fid_run)