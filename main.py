#!/usr/bin/env python
# -*- coding: latin-1 -*-
'''Function to evaluate the activity of each species during the simulation, 
   catalyst substrate product or nothing
'''

import sys, os # Standard library
import datetime as dt
import linecache as lc
from copy import deepcopy
import glob
from argparse import ArgumentParser
import numpy as np # Scientific library
from numpy import * 

try:
    from pylab import *
except:
    pass
   
from lib.IO import * 
from lib.graph import raf
from lib.graph import network
from lib.dyn import dynamics as dm

#ÊInput parameters definition 
if __name__ == '__main__':
	parser = ArgumentParser(
				description='Main script of ACS analysis.'
				, epilog='''ACS ANALYSIS Main File. ''') 
	parser.add_argument('-i', '--initanal', type=int, help='Analysis of the initial structures (def: 0)', choices=[0,1], default=0)
	parser.add_argument('-e', '--exhaustive', type=int, help='Analysis of the dynamic structures (def: 1)', choices=[0,1], default=1)
	parser.add_argument('-d', '--decay', type=float, help='Decay time (def: 0)', default=0)
	parser.add_argument('-t', '--timeWindow', type=float, help='Dynamical time window (def: 10 seconds)', default=10)
	parser.add_argument('-m', '--maxDim', help='Max Dimension of the system (def: 4)', default='4', type=int)
	parser.add_argument('-p', '--strPath', help='Path where files are stored (def: ./)', default='./')
	parser.add_argument('-r', '--resFolder', help='Name of the result folder (def: res)', default='res')
	parser.add_argument('-b', '--debug', type=int, help='Debug Mode (Def: 0)', default=0)
	
	args = parser.parse_args()
	
	# Create absolute paths
	strPath = os.path.abspath(args.strPath)
	tmpDirs = sort(os.listdir(strPath))
	
	# Goes into the simulation folder
	os.chdir(strPath)
	
	# Create stas folders
	ndn = '_0_new_allStatResults'
	newdirAllResults = os.path.join(strPath, ndn)
	if not os.path.isdir(newdirAllResults):
		try: os.mkdir(newdirAllResults)
		except: print "Impossible to create statistic directory", newdirAllResults; sys.exit(1)
	print "\n\n********************\n\n|- Simulation Folder: ", strPath
	# System Type
	_CLOSE_ = 0
	_PROTO_ = 1
	_CSTR_ = 2
	
	if args.initanal == 1:
		fname_initRafRes = os.path.join(newdirAllResults, '0_initRafAnalysis.csv')
		fname_initRafResLIST = os.path.join(newdirAllResults, '0_initRafAnalysisLIST.csv')
		fname_initRafResALL = os.path.join(newdirAllResults, '0_initRafAnalysisALL.csv')
		fid_initRafRes = open(fname_initRafRes, 'w')
		fid_initRafResLIST = open(fname_initRafResLIST, 'w')
		fid_initRafResALL = open(fname_initRafResALL, 'w')
		strToWrite = "Folder\tP\tAC\tM\tRAFsize\tClosureSize\tCatsSize\tuRAF\n"
		fid_initRafRes.write(strToWrite)
	
	for tmpDir in tmpDirs:
		
		totDirName = os.path.join(strPath,tmpDir)
		if os.path.isdir(totDirName):
			# Move to the directory 
			os.chdir(totDirName)
			
			# SIMULATION DYNAMICS FOLDER
			resDirPath = os.path.abspath(os.path.join("./", args.resFolder))
			print " \- Results Folder: ", totDirName
			if os.path.isdir(resDirPath):
				
				# Analysis of the initial structures 
				conf = readfiles.readConfFile(totDirName) # Configuration file upload
				# System type recognition
				if (conf[6] == 0) & (conf[7] > 0): sysType = _PROTO_
				elif (conf[6] > 0) & (conf[7] == 0): sysType = _CSTR_
				elif (conf[6] == 0) & (conf[7] == 0): sysType = _CLOSE_
				foodList = dm.generateFluxList(totDirName, sysType)
								
				# Initial Analysis are turned ON (RAF ANALYSIS)
				if args.initanal == 1:
					print "   |- LOADING init structures..."
					rcts = readfiles.loadAllData(totDirName,'_acsreactions.csv') # reaction file upload
					cats = readfiles.loadAllData(totDirName,'_acscatalysis.csv') # catalysis file upload
					# REAL RAF COMPUTATION (!!!!!! 1 is average connectivity)
					raf.rafComputation(fid_initRafRes, fid_initRafResALL, fid_initRafResLIST, tmpDir, conf[9], 1, rcts, cats, foodList, args.maxDim, debug=args.debug)
					
				os.chdir(resDirPath)
				print "    \-Oucomes Folder ", resDirPath
				# Find the number of generations
				numberOfGen = len(glob.glob(os.path.join(resDirPath,'times_*')))
				# For each generation
				os.chdir(resDirPath)
				
				# Analysis of the dynamics
				if (args.exhaustive == 1) | (args.decay > 0):
					print "\t|- Analysis of the dynamics..."
					for ngen in range(1,numberOfGen+1):
						print "\t|-Generation ", ngen
						strZeros = readfiles.zeroBeforeStrNum(ngen, numberOfGen)
						# IF CATALYIS 
						if args.exhaustive == 1:
							# Saving File Creation
							if args.exhaustive == 1:
								print "\t\t|- Analysis of the ongoing data structures... "
								fName = 'RAF_structuresInTime_analysis_gen_' + strZeros + str(ngen) + '.csv'
								fname_inTimeRafRes = os.path.join(newdirAllResults, fName)
								fid_inTimeRafRes = open(fname_inTimeRafRes, 'w')
								if conf[10] >= 1: 
									strToWrite = "t\t#CL\t#RAF\tRAF\n"
									potential = False
								else: 
									strToWrite = "t\t#CL\t#RAF\t#PCL\t#PRAF\tRAF\n"
									potential = True
								fid_inTimeRafRes.write(strToWrite)
							# reactions file loading
							if ngen == 1:
						  		strRctZero = 'reactions_' + strZeros + str(0) + '*'
						  		strCatZero = 'catalysis_' + strZeros + str(0) + '*'
						  		rctFilesZero = sorted(glob.glob(os.path.join(resDirPath,strRctZero)))
						  		catFilesZero = sorted(glob.glob(os.path.join(resDirPath,strCatZero)))
							strRct = 'reactions_' + strZeros + str(ngen) + '*'  
							strCat = 'catalysis_' + strZeros + str(ngen) + '*'  
						  	# Searching for files
						  	rctFiles = sorted(glob.glob(os.path.join(resDirPath,strRct)))
						  	catFiles = sorted(glob.glob(os.path.join(resDirPath,strCat)))
						  	
						  	if ngen == 1: # If the generation is the first one the 0000 files must be uploaded too
					  			rctFiles = rctFilesZero + rctFiles
					  			catFiles = catFilesZero + catFiles
					  		
					  		sngTime = conf[2] / (float(len(rctFiles)) - 1)
						  	actTime = 0
					  		
					  		for id, sngFile in enumerate(rctFiles):
					  			# Load files
					  			rcts = readfiles.loadAllData(totDirName,sngFile) # reaction file upload
								cats = readfiles.loadAllData(totDirName,catFiles[id]) # catalysis file upload
								foodList = dm.generateFluxList(totDirName, sysType)
								# Reshape rcts according to the reactions really occurred
								procrcts = rcts[rcts[:,5] > 0,:]
								proccats = cats[cats[:,3] > 0,:]
								R = raf.rafDynamicComputation(fid_inTimeRafRes, actTime, procrcts[:,0:5], proccats[:,0:5], foodList, potential, rcts, cats, debug=args.debug)
								#print R
								actTime += sngTime
						else:
							# Last reaction File is however loaded
							strRct = 'reactions_' + strZeros + str(ngen) + '*'  
							strCat = 'catalysis_' + strZeros + str(ngen) + '*'  
						  	# Searching for files
						  	rctFiles = sorted(glob.glob(os.path.join(resDirPath,strRct)))
						  	catFiles = sorted(glob.glob(os.path.join(resDirPath,strCat)))
						  	# Upload latest file
						  	lastRct = readfiles.loadAllData(totDirName,rctFiles[-1]) # reaction file upload
						  	lastCat = readfiles.loadAllData(totDirName,catFiles[-1]) # reaction file upload
							
						# If decay analysis is TURNED ON (rections are added step by step to the graph and they are removed 
						# if they do not occur twice within the decay time
						if args.decay > 0:
							print "\t\t|- Analysis of the reactions actually occurred according to the decay time equal to ", str(args.decay)
							fName = 'RAF_dynamics_analysis_gen_' + strZeros + str(ngen) + '.csv'
							fname_dynRafRes = os.path.join(newdirAllResults, fName)
							fid_dynRafRes = open(fname_dynRafRes, 'w')
							
							# FILE FIRST RAW
							strToWrite = "t\t#CL\t#RAF\tRAF\n"
							fid_dynRafRes.write(strToWrite)
							
							strRctPar = 'reactions_parameters_' + strZeros + str(ngen) + '*'
							rctParamFile = sorted(glob.glob(os.path.join(resDirPath,strRctPar)))
							fid = open(rctParamFile[0], 'r')
							previousTime = 0
							decayTime = args.decay
							condensation_counter = 0
							endo_condensation_counter = 0
							cleavage_counter = 0
							endo_cleavage_counter = 0
							nAnal = 1
							rctCurrID = 0
							catCurrID = 0
							for rctL, line in enumerate(fid):								 
								# Load single reaction parameters
								reaction, rtime, cc, cat, mol_I, mol_II, mol_III, loadedMolsConc, loadedMols,\
		   						gillMean, gillSD, gillEntropy, newSpeciesCreationProb, reverseProbability = readfiles.splitRctParsLine(line)
								
								if rctL % 10000 == 0: print "\t\t\t|- Reaction number ", rctL, " - Time: ", rtime
								
								# update time intervals 
								timeInterval = rtime - previousTime
								previousTime = rtime
				
								if rctL > 0:
									# REACTIONS NOT OCCURING TWICE IN THE DACAY TIME ARE REMOVED FROM THE GRAPH
									graph = network.removeRareRcts(graph,2,3,4,timeInterval)
									graphSUB = network.removeRareRcts(graphSUB,2,3,4,timeInterval)
									onrcts = network.removeRareRcts(onrcts,5,6,7,timeInterval)
									oncats = network.removeRareRcts(oncats,3,4,5,timeInterval)
								
								# ONGOING REACTION STRUCTURES CREATION
								if rctL == 0:
									mol_I, mol_II, mol_III = network.fixCondensationReaction(mol_I, mol_II, mol_III, lastRct)
									onrcts = np.array([[rctCurrID, cc, mol_I, mol_II, mol_III, decayTime, 0, decayTime, 1]], dtype=np.float64)
									oncats = np.array([[catCurrID, cat, onrcts[0,0], decayTime, 0, decayTime, 1]], dtype=np.float64)
									rctCurrID += 1
									catCurrID += 1
								else:
									# !!! ONCSTR ANALYSIS
									mol_I, mol_II, mol_III = network.fixCondensationReaction(mol_I, mol_II, mol_III, lastRct)
									if sum((onrcts[:,1] == cc) & (onrcts[:,2] == mol_I) & (onrcts[:,3] == mol_II)) == 1:
										positionR = ((onrcts[:,1] == cc) & (onrcts[:,2] == mol_I) & (onrcts[:,3] == mol_II))	
										onrcts[positionR] = [onrcts[positionR,0], cc, mol_I, mol_II, mol_III, decayTime, 0, decayTime, onrcts[positionR,8]+1]	
									else:
		
										onrcts = np.vstack([onrcts,(rctCurrID, cc, mol_I, mol_II, mol_III, decayTime, 0, decayTime, 1)])	
										positionR = onrcts[:,0] == rctCurrID
										rctCurrID += 1	
																		
									
									# !!! ONCAT ANALYSIS
									if sum((oncats[:,1] == cat) & (oncats[:,2] == onrcts[positionR,0])) == 1:
										#print  onrcts[positionR,0]
										position = ((oncats[:,1] == cat) & (oncats[:,2] == onrcts[positionR,0]))	
										oncats[position] = [oncats[position,0], cat, onrcts[positionR,0], decayTime, 0, decayTime, oncats[position,6]+1]	
									else:
										oncats = np.vstack([oncats,(catCurrID, cat, onrcts[positionR,0], decayTime, 0, decayTime, 1)])
										catCurrID += 1	
								
								# RAF ANALYSIS	
								if rtime >= float((args.timeWindow * nAnal)):
									print "\t\t\t|- RAF analysis...", " ", rtime, " - Reaction ", rctL
									foodList = dm.generateFluxList(totDirName, sysType)
									R = raf.rafDynamicComputation(fid_dynRafRes, rtime, onrcts[:,0:5], oncats[:,0:5], foodList, False, completeRCTS=lastRct , debug=args.debug)
									nAnal += 1
								
								# REACTION GRAPH CREATION
								
								if (cc == 0)|(cc == 7): # CONDENSATION and ENDO CONDENSATION
									if cc == 0:
										condensation_counter += 1
									else:
										endo_condensation_counter +=  1
										
									if rctL == 0: #If it is the first reaction nparray (matrix) is created
										graph = np.array([[cat, mol_I, decayTime, float(0), decayTime, 1]])
										graphSUB = np.array([[mol_II, mol_I, decayTime, float(0), decayTime, 1]])
										if mol_II != mol_III:
											graphSUB = np.vstack([graphSUB,(mol_III, mol_I, decayTime, 0, decayTime, 1)]) # Substrate 2 (If different from 1)
									else: 
										#!!! CAT -> PRO, Otherwise if the reaction is already present its parameters are updated
										if sum((graph[:,0] == cat) & (graph[:,1] == mol_I)) == 1:
											position = ((graph[:,0] == cat) & (graph[:,1] == mol_I))
											graph[position] = [cat, mol_I, decayTime, 0, decayTime, graph[position,5]+1]
										else:
											# Otherwise a new reaction is added at the end of the matrix
											graph = np.vstack([graph,(cat, mol_I, decayTime, 0, decayTime, 1)])
							
										#!!! SUB -> PRO, Otherwise if the reaction is already present its parameters are updated
										if sum((graphSUB[:,0] == mol_II) & (graphSUB[:,1] == mol_I)) == 1:
											position = ((graphSUB[:,0] == mol_II) & (graphSUB[:,1] == mol_I))
											graphSUB[position] = [mol_II, mol_I, decayTime, 0, decayTime, graphSUB[position,5]+1]
										else:
											# Otherwise a new reaction is added at the end of the matrix
											graphSUB = np.vstack([graphSUB,(mol_II, mol_I, decayTime, 0, decayTime, 1)])
										if mol_II != mol_III:
											if sum((graphSUB[:,0] == mol_III) & (graphSUB[:,1] == mol_I)) == 1:
												position = ((graphSUB[:,0] == mol_III) & (graphSUB[:,1] == mol_I))
												graphSUB[position] = [mol_III, mol_I, decayTime, 0, decayTime, graphSUB[position,5]+1]
											else:
												graphSUB = np.vstack([graphSUB,(mol_III, mol_I, decayTime, 0, decayTime, 1)])	
														
										
								else: # CLEAVAGE (The same of condensation expected for the double computation due to the presence of two products
									if cc == 6:
										endo_cleavage_counter += 1
									else:
										cleavage_counter +=  1
										
									if rctL == 0: # If it is the first reaction nparray is created
										# CAT -> PROD
										graph = np.array([[cat, mol_II, decayTime, 0, decayTime, 1]]) # Product 1
										if mol_II != mol_III:
											graph = np.vstack([graph,(cat, mol_III, decayTime, 0, decayTime, 1)]) # Product 2 (If different from 1)
										# SUB -> PRO
										graphSUB = np.array([[mol_I, mol_II, decayTime, 0, decayTime, 1]])
										if mol_II != mol_III:
											graphSUB = np.vstack([graphSUB,(mol_I, mol_III, decayTime, 0, decayTime, 1)]) # Sub 2 (If different from 1)					
									else:
										if sum((graph[:,0] == cat) & (graph[:,1] == mol_II)) == 1:
											position = ((graph[:,0] == cat) & (graph[:,1] == mol_II))
											graph[position] = [cat, mol_II, decayTime, 0, decayTime, graph[position,5]+1]
										else:
											graph = np.vstack([graph,(cat, mol_II, decayTime, 0, decayTime, 1)])				
										if mol_II != mol_III:
											if sum((graph[:,0] == cat) & (graph[:,1] == mol_III)) == 1:
												position = ((graph[:,0] == cat) & (graph[:,1] == mol_III))
												graph[position] = [cat, mol_III, decayTime, 0, decayTime, graph[position,5]+1]
											else:
												graph = np.vstack([graph,(cat, mol_III, decayTime, 0, decayTime, 1)])
										# SUB -> PRO, Otherwise if the reaction is already present its parameters are updated
										if sum((graphSUB[:,0] == mol_I) & (graphSUB[:,1] == mol_II)) == 1:
											position = ((graphSUB[:,0] == mol_I) & (graphSUB[:,1] == mol_II))
											graphSUB[position] = [mol_I, mol_II, decayTime, 0, decayTime, graphSUB[position,5]+1]
										else:
											# Otherwise a new reaction is added at the end of the matrix
											graphSUB = np.vstack([graphSUB,(mol_I, mol_II, decayTime, 0, decayTime, 1)])
										if mol_II != mol_III:
											if sum((graphSUB[:,0] == mol_I) & (graphSUB[:,1] == mol_III)) == 1:
												position = ((graphSUB[:,0] == mol_I) & (graphSUB[:,1] == mol_III))
												graphSUB[position] = [mol_I, mol_III, decayTime, 0, decayTime, graphSUB[position,5]+1]
											else:
												graphSUB = np.vstack([graphSUB,(mol_I, mol_III, decayTime, 0, decayTime, 1)])			
								
							
							
								#raw_input("...")
							fid.close() # close fid
					
					
	if args.initanal == 1: 
		fid_initRafRes.close()
		fid_initRafResLIST.close()
		fid_initRafResALL.close()
	
	if args.exhaustive == 1: fid_inTimeRafRes.close()
	
					
				  	
				  	


