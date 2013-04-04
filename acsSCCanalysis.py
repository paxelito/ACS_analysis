#!/usr/bin/env python
# -*- coding: latin-1 -*-

'''Function to analyze the different attractors emerging from different simulations in time. 
'''

import sys, os # Standard library
import datetime as dt
import glob
import numpy as np # Scientific library
from numpy import * 
import networkx as nx

try:
    from pylab import *
except:
    pass

# Function to create string zero string vector before graph filename. 
# According to the total number of reactions N zeros will be add before the instant reaction number 
# (e.g. reaction 130 of 10000 the string became '00130')
def zeroBeforeStrNum(tmpl, tmpL):
	strZero = ''
	nZeros = len(str(tmpL)) - len(str(tmpl))
	if nZeros > 0:
		for i in range(0,nZeros): strZero = strZero + '0'
	return strZero
#eof
	
#Function to load the reaction graph
def loadReactionGraph():
	os.chdir(StrPath)
	os.chdir(os.path.join(StrPath,folderCat))
	outFname = '_iGraph_CAT' + filextPre + "*"
	print os.path.join(StrPath,folderCat,outFname)
	files = sorted(glob.glob(outFname))
	print files
	try:
		fid = open(files[len(files)-1], 'r')
	except:
		print 'impossible to load ', files; sys.exit(1)	
		
	print "\tloading reactions graph..."
	rctID = 1
	for line in fid:
		tmpCat, tmpProd, tmpdt1, tmpdt2, tmpdt3, tmpTimes = line.split()
		cat = int(tmpCat)
		prod = int(tmpProd)
		dt1 = float(tmpdt1)
		dt2 = float(tmpdt2)
		dt3 = float(tmpdt3)
		t = int(tmpTimes)
		if rctID == 1:
			graph = np.array([[cat, prod, dt1, dt2, dt3, t]])
			rctID = 0
		else:
			graph = np.vstack([graph,(cat, prod, dt1, dt2, dt3, t)])
			
	os.chdir(StrPath)
	
	return graph
	
#Function to load a specific reaction graph (defined as a parameter)
def loadSpecificReactionGraph():
	os.chdir(StrPath)
	print "Load reactions file: ", tmpRctFileToLoad
	try:
		fid = open(tmpRctFileToLoad, 'r')
	except:
		print 'impossible to load ', files; sys.exit(1)	
		
	print "\tloading reactions graph..."
	rctID = 1
	for line in fid:
		tmpCat, tmpProd, tmpdt1, tmpdt2, tmpdt3, tmpTimes = line.split()
		cat = int(tmpCat)
		prod = int(tmpProd)
		dt1 = float(tmpdt1)
		dt2 = float(tmpdt2)
		dt3 = float(tmpdt3)
		t = int(tmpTimes)
		if rctID == 1:
			graph = np.array([[cat, prod, dt1, dt2, dt3, t]])
			rctID = 0
		else:
			graph = np.vstack([graph,(cat, prod, dt1, dt2, dt3, t)])
						
	os.chdir(StrPath)
	
	return graph
	
#Function to load a specific reaction SUBgraph (defined as a parameter)
def loadSpecificReactionSubGraph():
	os.chdir(StrPath)
	print "Load reactions file: ", tmpRctSUBFileToLoad
	try:
		fid = open(tmpRctSUBFileToLoad, 'r')
	except:
		print 'impossible to load ', files; sys.exit(1)	
		
	print "\tloading reactions graph..."
	rctID = 1
	for line in fid:
		tmpSub, tmpProd, tmpdt1, tmpdt2, tmpdt3, tmpTimes = line.split()
		sub = int(tmpSub)
		prod = int(tmpProd)
		dt1 = float(tmpdt1)
		dt2 = float(tmpdt2)
		dt3 = float(tmpdt3)
		t = int(tmpTimes)
		if rctID == 1:
			graphSUB = np.array([[sub, prod, dt1, dt2, dt3, t]])
			rctID = 0
		else:
			graphSUB = np.vstack([graphSUB,(sub, prod, dt1, dt2, dt3, t)])
						
	os.chdir(StrPath)
	
	return graphSUB
	
#Function to save graphs to file (Considering that all variables all global in reading)
def saveGraphToFile():
	os.chdir(StrPath)
	tmpStrZeros = zeroBeforeStrNum(reaction, totalRcts)
	outFname = '_iGraph_CAT' + filextPre + '_' + tmpStrZeros + str(reaction) + '_' + str(rtime) + '.csv'
	saveFile = open(outFname, 'w')
	# Move to the directory 
	np.savetxt(saveFile, graph, '%d\t%d\t%10.5f\t%10.5f\t%10.5f\t%d')
	saveFile.close()
	os.rename(outFname, os.path.join(folderCat,outFname))

#Function to save graphs SUB->PRO to file (Considering that all variables all global in reading)
def saveGraphSUBToFile():
	os.chdir(StrPath)
	tmpStrZeros = zeroBeforeStrNum(reaction, totalRcts)
	outFname = '_iGraph_SUB' + filextPre + '_' + tmpStrZeros + str(reaction) + '_' + str(rtime) + '.csv'
	saveFile = open(outFname, 'w')
	# Move to the directory 
	np.savetxt(saveFile, graphSUB, '%d\t%d\t%10.5f\t%10.5f\t%10.5f\t%d')
	saveFile.close()
	os.rename(outFname, os.path.join(folderSub,outFname))

#Function to save energy file 
def saveNrgToFile():
	os.chdir(StrPath)
	outFname = '_ENERGY.csv'
	saveFile = open(outFname, 'w')
	# Move to the directory 
	np.savetxt(saveFile, nrgTimeSeries, '%10.5f\t%10.5f\t%d')
	saveFile.close()
	os.rename(outFname, os.path.join(tmpDir,outFname))

#Function to save gillespie information files
def saveGillToFile():
	os.chdir(StrPath)
	outFname = '_GILLESPIE.csv'
	saveFile = open(outFname, 'w')
	# Move to the directory 
	np.savetxt(saveFile, gillTimeSeries, '%10.5f\t%10.5f\t%10.5f\t%10.5f')
	saveFile.close()
	os.rename(outFname, os.path.join(tmpDir,outFname))

# get path for placing simulation
# INPUT PARAMETERS
try:
	StrPath = sys.argv[1] # Here the path of the simulation output file
	threshold = float(sys.argv[2]) # Theshold above which the concentration has to be 
	decayTime = float(sys.argv[3]) # Reaction decay time
	tmpRctFileToLoad = sys.argv[4] # Reaction file to load
	tmpRctSUBFileToLoad = sys.argv[5] # Reaction SUB file to load
except:
	print "Usage:",sys.argv[0], "infile outfile"; sys.exit(1)
	
print "Simulation Results Path: ", StrPath

# Create absolute path 
StrPath = os.path.abspath(StrPath)
today = dt.date.today()
mswindows = (sys.platform == "win32")

# Create summary file
cmdFileName = StrPath + '/' + str(today) + '_' + str(threshold) + '_convOverThreshold_' + str(decayTime) + '.csv'
if  os.path.isfile(cmdFileName):
	os.remove(cmdFileName)
cmdFileFid = open(cmdFileName, 'a')
strToWrite = "\tReaction Probability" + "\tEnergy Concentration" + "\tMean Over Threshold" + "\tTot over Threshold" + "\tReal ACSs" + "\tin ACSs" + "\tLeaves" + "\tby Chain" + "\tOverlap" + "\tAutocatalysis" + "\tin ACSs (W)" + "\tLeaves (W)" + "\tby Chain (W)" + "\tOverlap (W)" + "\tAutocatalysis (W)" +"\t[ACS]" + "\t[Leaves]" + "\t[Chains]" + "\t[Self]" + "\tEndo Cond" + "\tCond" + "\tEndo Cleavage" + "\tCleavage" + "\tWasted species" + "\n"
cmdFileFid.write(strToWrite)

initilizeGraphStructure = 0;

# If the tmpRctFileToLoad is not empty before to begin the analysis the initial reaction file is loaded
if tmpRctFileToLoad != "":
	graph = loadSpecificReactionGraph()
	graphSUB = loadSpecificReactionSubGraph()
	initilizeGraphStructure = 1;
	tmpRctFileToLoadSplitString = tmpRctFileToLoad.split("_")
	# Retry the reaction number 
	initRctId = int(tmpRctFileToLoadSplitString[len(tmpRctFileToLoadSplitString)-2])
	initTempTime = tmpRctFileToLoadSplitString[len(tmpRctFileToLoadSplitString)-1]
	initTime = float(initTempTime[0:len(initTempTime)-4])
	print initRctId
	print initTime
else:
	initRctId = 1;
	initTime = 0;
	
	
# Create directories list 
tmpDirs = sort(os.listdir(StrPath))
	
# For each folder
for tmpDir in tmpDirs:
	print "*********************************************************************************"
	print "* NEW SIMULATION ANALYSIS PROCESS"
	print "* Folder: ",tmpDir
	print "*********************************************************************************"
	os.chdir(StrPath)
	totDirName = os.path.join(StrPath,tmpDir)
	if os.path.isdir(totDirName):
		# Move to the directory 
		os.chdir(totDirName)
		# if there is a result folder (i.e. totDirName is a simulation folder) 
		resDirPath = os.path.abspath(os.path.join("./", "res"))
		print "  |- Results Folder: ", resDirPath
		if os.path.isdir(resDirPath):
			
			os.chdir(resDirPath)
			
			# Find the number of generations contained in this simulation folder
			numberOfGen = len(glob.glob(os.path.join(resDirPath,'times_*')))
			
			# FOR EACH GENERATION
			for nGen in range(1,numberOfGen+1):
				
				os.chdir(resDirPath)
				
				# select zero chars to insert before the generation number according to the folder structure.
				strZeros = zeroBeforeStrNum(nGen, numberOfGen)
				# Open Parameter File
				paramFile = "acsm2s.conf"
				simFolder = 0
				try:
					fid = open(paramFile, 'r')
				except:
					print 'impossible to load ', paramFile; sys.exit(1)
								
				# Read reaction probability from configuration file
				for line in fid:
					strLine = line.split('=')
					if strLine[0] == "reactionProbability":
						rp = float(strLine[1])
					if strLine[0] == "nSeconds":
						if decayTime == 0:
							decayTime = int(strLine[1])
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
						maxLOut = float(strLine[1])
				
				# analysisTimeInterval is 1/10 of the total decay time
				_ANALTIMEINTERAVAL_ = totTimes / 10
				_ANALTIMEINTERAVALNOSAVE_ = totTimes / 100
						
				if nrgType < 2:
					nrg = 1
				else:
					nrg = 0
						
				# Influx in Protocell 
				if (influx_rate == 0) and (maxLOut > 0):
					print '   |- Createing influx for protocell simulation...'
					speciesInFlux = range(0,int(pow(2,(maxLOut+1)) - 2))
				else: # influx in CSTR
					try:
						# Create influx species list
						print '   |- Createing influx for CSTR simulation...'
						fidflux = open('_acsinflux.csv', 'r')
						speciesInFlux = []
						for fluxID in fidflux:
							tmpID, tmpProb = fluxID.split()
							speciesInFlux.append(int(tmpID))
						fidflux.close()
					except:
						print 'impossible to load ', influxFile;
						print 'Influx will be set to empty'
						speciesInFlux = []
				
				print "    <> Influx IDsm -> ", speciesInFlux	
				
				# initialize stat variables 
				overThreshold = float(0)
				overThresholdTOT = float(0)
				# Searching for the Species files related to the nGen-th generation
				tmpFilesToSearch = 'species_' + strZeros + str(nGen) + '_*'
				speciesFiles = sorted(glob.glob(os.path.join(resDirPath,tmpFilesToSearch)))
				speciesFile = speciesFiles[-1]
				
				print "   |- Species File: ", speciesFile, "--------"
				
				# Open Species File
				try:
					fidSpecies = open(speciesFile, 'r')
				except:
					print '   |- !!!impossible to load ', speciesFile; sys.exit(1)
				
				# For each last species file
				ok = 0
				IDsOverThreshold = []
				concVec = []
				
				# for each species compute species over threshold
				for sp in fidSpecies:
				
					tmpID, tmpSeq, tmpConc, tmpDiff, tmpSol, tmpCpxDiss, tmpCpxCut, tmpEval, tmpAge,\
					tmpReb, tmpCatID, tmpSubID, tmpKpho, tmpLoadConc, tmpConcFixed = sp.split()
					
					index = int(tmpID)
					conc = float(tmpConc)
					cpxCut = int(tmpCpxCut)
					age = float(tmpAge)
						
					#if the species is not in the influx and greater than the max length
					if index not in speciesInFlux:
						realT = threshold # To complete to be consistent
						concVec.append(conc)
						if (conc > realT) & (cpxCut == 0) & (age > 0):
							#print tmpSeq, " ", conc, int(len(tmpSeq))
							if ok == 0:
								overThreshold += 1
								ok = 1
							overThresholdTOT += 1
							IDsOverThreshold.append(int(tmpID))
					else:
						concVec.append(conc)
	
								
				# If there are species over threshold network analysis is performed	
				print "    |- Graph Creation"
				folderCat = '__0_iGraph_CAT_' + str(decayTime)
				folderSub = '__0_iGraph_SUB_' + str(decayTime)
				
				# Move to the root directory directory 
				os.chdir(StrPath)
				
				# Create iGraph directory where iGraph_CAT file will be stored
				# If the folder has been already create, it is removed
				newdir = os.path.join(os.curdir, folderCat)
				if not os.path.isdir(newdir):
					try:
						os.mkdir(newdir)
					except:
						print "   |- !!! Impossible to create iGraph directory"; sys.exit(1)
			
				# Create iGraph directory where iGraph_SUB file will be stored. 
				# If the folder has been already create, it is removed
				newdirSUB = os.path.join(os.curdir, folderSub)
				if not os.path.isdir(newdirSUB):
					try:
						os.mkdir(newdirSUB)
					except:
						print "   |- !!! Impossible to create iGraph directory"; sys.exit(1)
				
				if os.path.isdir(os.path.join(tmpDir, "res")):
					os.chdir(os.path.join(tmpDir, "res"))		
					
					
					filextPre = '_' + zeroBeforeStrNum(nGen, numberOfGen) + str(nGen)
			
					# Select the reaction_parameters file according to the generation
					rctParamFileQ = 'reactions_parameters_' + strZeros + str(nGen) + '_*'
					rctParamFile = sorted(glob.glob(os.path.join(resDirPath,rctParamFileQ)))
					#open file
					fid = open(rctParamFile[0], 'r')
					rctIDshow = 1.0
					rctIDshowNoSave = 1.0
					rctID = initRctId
					previousTime = 0
					
					#Initialize reaction type counter
					endo_condensation_counter = 0;
					condensation_counter = 0;
					endo_cleavage_counter = 0;
					cleavage_counter = 0;
					
					print "   |- File ", rctParamFile[0], " processing..."
										
					if nrg == 1:
						nrgTimeSeries = np.array([[0,0,0]])
						
					gillTimeSeries = np.array([[0,0,0,0]])
					
					tmpFlagLastRctSaved = 0 # Flag var to save the last reaction
					for line in fid:
											
						# Load single reaction parameters
						tmpReaction, tmpTime, tmpcc, tmpCat, tmpMol_I, tmpMol_II, tmpMol_III, tmpLoadedMols,\
						 tmpLoadedMolsConc, tmpGillMean, tmpGillSD, tmpGillEntropy, tmpNSCprob, tmpRevProb = line.split()
						 
						reaction = int(tmpReaction) + initRctId
						rtime = float(tmpTime) + initTime
						cc = int(tmpcc) 
						cat = int(tmpCat)
						mol_I = int(tmpMol_I)
						mol_II = int(tmpMol_II)
						mol_III = int(tmpMol_III)
						loadedMolsConc = float(tmpLoadedMolsConc)
						loadedMols = int(tmpLoadedMols)
						gillMean = float(tmpGillMean)
						gillSD = float(tmpGillSD)
						gillEntropy = float(tmpGillEntropy)
						newSpeciesCreationProb = float(tmpNSCprob)
						reverseProbability = float(tmpRevProb)
		
						if nrg == 1:
							nrgTimeSeries = np.vstack([nrgTimeSeries,(rtime, loadedMolsConc, loadedMols)])
							
						gillTimeSeries = np.vstack([gillTimeSeries,(rtime, gillMean, gillSD, gillEntropy)])
													
						printTemporalMessage = 1
						# If the time is right save iGraph structures to file
						if (rtime >= rctIDshow * _ANALTIMEINTERAVAL_):
							print "- ", tmpDir, " | Gen", nGen, " | Reaction ", reaction, " - time: ", rtime, " - Save structures"
							saveGraphToFile()
							saveGraphSUBToFile()
							rctIDshow = rctIDshow + 1
							printTemporalMessage = 0
							
						# If the system is close to the end graphs are saved.
# 						if (tmpFlagLastRctSaved == 0) and (rtime > (totTimes - 1)) and (rtime < totTimes):
# 							print "- ", tmpDir, " | Gen", nGen, " | Reaction ", reaction, " - time: ", totTimes, " - Save structures"
# 							saveGraphToFile()
# 							saveGraphSUBToFile()
# 							tmpFlagLastRctSaved = 1
							 
			
						# If the prompt time is right a message on the screen indicating the reaction number and the time is shown
						if ((rtime > rctIDshowNoSave * _ANALTIMEINTERAVALNOSAVE_) and (printTemporalMessage == 1)):
							print "- ", tmpDir, " | Gen", nGen, " | Reaction ", reaction, " - time: ", rtime
							rctIDshowNoSave = rctIDshowNoSave + 1
			
						# update time intervals 
						timeInterval = rtime - previousTime
						previousTime = rtime
		
						# Increment the time of each reaction present in the system
						# If the time of some reactions overcome the decay time the reaction is removed from the igraph structure
						if rctID > 1:
							# CATALYST -> PRODUCT
							if shape(graph)[0] > 1:
								graph[:,3] = graph[:,3] + timeInterval
								graph[:,4] = graph[:,2] - graph[:,3]
								graph = graph[graph[:,4]>0,:]
							# SUBSTRATE -> PRODUCT
							if shape(graphSUB)[0] > 1:
								graphSUB[:,3] = graphSUB[:,3] + timeInterval
								graphSUB[:,4] = graphSUB[:,2] - graphSUB[:,3]
								graphSUB = graphSUB[graphSUB[:,4]>0,:]
				
						if (cc == 0)|(cc == 7): # CONDENSATION and ENDO CONDENSATION
							if cc == 0:
								condensation_counter += 1
							else:
								endo_condensation_counter +=  1
								
							if rctID == 1: #If it is the first reaction nparray (matrix) is created
								graph = np.array([[cat, mol_I, decayTime, 0, decayTime, 1]])
								graphSUB = np.array([[mol_II, mol_I, decayTime, 0, decayTime, 1]])
								if mol_II != mol_III:
									graphSUB = np.vstack([graphSUB,(mol_III, mol_I, decayTime, 0, decayTime, 1)]) # Substrate 2 (If different from 1)
							else: 
								# CAT -> PRO, Otherwise if the reaction is already present its parameters are updated
								if sum((graph[:,0] == cat) & (graph[:,1] == mol_I)) == 1:
									position = ((graph[:,0] == cat) & (graph[:,1] == mol_I))
									graph[position] = [cat, mol_I, decayTime, 0, decayTime, graph[position,5]+1]
								else:
									# Otherwise a new reaction is added at the end of the matrix
									graph = np.vstack([graph,(cat, mol_I, decayTime, 0, decayTime, 1)])
					
								# SUB -> PRO, Otherwise if the reaction is already present its parameters are updated
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
								
							if rctID == 1: # If it is the first reaction nparray is created
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
						
						# Reactions counters is updated			
						rctID = rctID + 1	
					# Save Gillespie structure on file 	
					saveGillToFile()
						
					print "\n |-  ------- NetworkX analysis ---------"
								
					# Network analysis
					# Initialize graph

					print "   |- digraph creation... "
					realSccs = 0
					Gcatpro = nx.DiGraph()
					for sngRct in graph:
						Gcatpro.add_weighted_edges_from([(sngRct[0],sngRct[1],sngRct[5])])

					print "   |- Strongly connected components analysis... "
					
					scc = nx.strongly_connected_components(Gcatpro)
					sccN= nx.number_strongly_connected_components(Gcatpro)
					selfLoops = Gcatpro.number_of_selfloops()
					selfLoopsEgdes = Gcatpro.selfloop_edges()
					realSccs += selfLoops
					
					# Computing actual number od strongly connected components
					prod_inSCC = 0
					prod_chain = 0
					prod_bySCC = 0
					prod_overlap = 0
					sccID = 0
					autocatalysis = 0
					prod_inSCC_weight = 0
					prod_chain_weight = 0
					prod_bySCC_weight = 0
					prod_overlap_weight = 0
					self_loop_weight = 0
					conc_inSCC = 0
					conc_chain = 0
					conc_bySCC = 0
					conc_overLap = 0
					conc_selfCat = 0
					wasteSpecies = 0
					
					for c in scc: 
						if len(c) > 1:
							sccID += 1
							print "<> SCC n.", sccID
							print c
							realSccs += 1

					print "   |- Number of ACS:", realSccs - selfLoops
					print "   |- Number of ACS(length 1):", selfLoops, "..."
					
					print "    <> IDs over threshold -> ", IDsOverThreshold
					
					if len(IDsOverThreshold) > 0: # If there are species over threshold SOT
						for IdsOT in IDsOverThreshold: # For each species over threshold
							alreadyAdded_ACS = 0
							alreadyAdded_leaves = 0
							alreadyAdded_chain = 0
							tmpProd_chain = 0
							try: # If there are incoming nodes 
								incomingNodes = Gcatpro.predecessors(IdsOT) # incoming nodes
								tmpProd_chain += len(incomingNodes) # Number of incoming nodes 
								tempProd_chain_weight = 0
								if len(incomingNodes) > 0:
									for sngInNode in incomingNodes: # for each single incoming nodes Xj -> Xi
										alreadyAdded_ACS = 0
										alreadyAdded_leaves = 0
										noInAcs = 1 
										if realSccs > 0: # If there are ACS
											for c in scc:  # for each ACS (if present)
												if len(c) > 1:
													inSCCFlag = 0 
													if IdsOT in c: # If the SOT is in an ACS
														inSCCFlag = 1
														if alreadyAdded_ACS == 0:
															conc_inSCC += concVec[IdsOT]
															alreadyAdded_ACS = 1
													weightToDistribute = graph[((graph[:,0] == sngInNode) & (graph[:,1] == IdsOT)),5]
													if sngInNode in c: # If the incoming node is in the ACS
														noInAcs = 0
														if inSCCFlag == 1: # If the node is in an ACS, then it has been produced within an ACS
															prod_inSCC += 1
															tmpProd_chain -= 1
															print "\t Within Acs", sngInNode, "->", IdsOT, "#", weightToDistribute, "-[IdsOT]",IdsOT, concVec[IdsOT]
															prod_inSCC_weight += weightToDistribute
														else: # Otherwise it has been produced by an ACS, so it is a leaves of the first degree
															if alreadyAdded_leaves == 0:
																conc_bySCC += concVec[IdsOT]
																alreadyAdded_leaves = 1
															prod_bySCC += 1
															tmpProd_chain -= 1
															print "\t From Acs", sngInNode, "->", IdsOT, "#", weightToDistribute, "-[IdsOT]",IdsOT, concVec[IdsOT]
															prod_bySCC_weight += weightToDistribute
	
										if noInAcs == 1: # If the species is not in the ACS and the incoming node as well
											prod_chain_weight += graph[((graph[:,0] == sngInNode) & (graph[:,1] == IdsOT)),5]
											if alreadyAdded_chain == 0:
												conc_chain += concVec[IdsOT]
												alreadyAdded_chain = 1
												
								if selfLoops > 0:
									for sngSelfLoop in selfLoopsEgdes:
										if IdsOT in sngSelfLoop:
											self_loop_weight +=  graph[((graph[:,0] == IdsOT) & (graph[:,1] == IdsOT)),5]
											conc_selfCat += concVec[IdsOT]
											autocatalysis += 1
							except:
								wasteSpecies += 1 # The specie is still over the threshold but no longer produced by the system into the considered time window
										
														
							prod_overlap += min(prod_bySCC, prod_inSCC, tmpProd_chain)	
							prod_chain += tmpProd_chain
							prod_overlap_weight += min(prod_inSCC_weight, prod_bySCC_weight, prod_chain_weight)
									
									
					print "   |- <> Number of Structural Autocatalytic set of molecules:", realSccs
					print "   |- <> -----------------------------------------------"
					print "   |- <> Species over threshold produced by a CHAIN:", prod_chain
					if realSccs > 0:
						print "   |- <> Species over threshold produced INTO an ACS:", prod_inSCC
						print "   |- <> Species over threshold produced INTO an ACS (weigthed):", prod_inSCC_weight
						print "   |- <> Species over threshold produced BY an ACS:", prod_bySCC
						print "   |- <> Species over threshold produced BY an ACS (weigthed):", prod_bySCC_weight
						print "   |- <> Species over threshold produced by an overlap:", prod_overlap
						print "   |- <> Species over threshold produced by an overlap (weigthed):", prod_overlap_weight
						print "   |- <> Species over threshold produced by itself:", autocatalysis
						print "   |- <> Species over threshold produced by itself (weigthed):", self_loop_weight
						print "   |- <> Concentration in ACSs:", conc_inSCC
						print "   |- <> Concentration in ACSs leaves:", conc_bySCC
						print "   |- <> Concentration in chains:", conc_chain
						print "   |- <> Concentration of autocatalyst:", conc_selfCat
						print "   |- <> Number of endo condensations:", endo_condensation_counter
						print "   |- <> Number of condensations:", condensation_counter
						print "   |- <> Number of endo cleavages:", endo_cleavage_counter
						print "   |- <> Number of cleavages:", cleavage_counter
					print "      <> -----------------------------------------------"						
							
					inDegreeMean = mean(Gcatpro.in_degree().values())
							
					print "   |- <> Mean IN-DEGREE:", inDegreeMean
			
					#print graph at the end of the analysis (the final situation is always stored)
					print "   |- <> File ", rctParamFile[0], " has been processed. Final Dimension: ", graph.shape[0]
					
					fid.close()
					saveGraphToFile()
					saveGraphSUBToFile()
							
			     				
					
				# Close species file
				fidSpecies.close()
					
				print "   |- <> ID over threshold (influx is not taken into account)"
				print IDsOverThreshold
				print overThreshold
				print overThresholdTOT
				print "   |- Species file length", float(len(speciesFiles))
				meanOverThreshold = float(overThreshold) / float(len(speciesFiles))
											
				strToWrite = tmpDir + "\t" + str(rp) + "\t" + str(nrgConc) + "\t" +  str(meanOverThreshold) + "\t" + str(overThresholdTOT) + "\t" 
				strToWrite += str(realSccs) + "\t" + str(prod_inSCC) + "\t" + str(prod_bySCC) + "\t" + str(prod_chain) + "\t" + str(prod_overlap) + "\t" + str(autocatalysis) + "\t"
				strToWrite += str(int(prod_inSCC_weight)) + "\t" + str(int(prod_bySCC_weight)) + "\t" + str(int(prod_chain_weight)) + "\t" + str(int(prod_overlap_weight)) + "\t" + str(int(self_loop_weight)) +"\t"
				strToWrite += str(conc_inSCC) + "\t" + str(conc_bySCC) + "\t" + str(conc_chain) + "\t" + str(conc_selfCat) + "\t"
				strToWrite += str(endo_condensation_counter) + "\t" + str(condensation_counter) + "\t" + str(endo_cleavage_counter) + "\t" + str(cleavage_counter) + "\t" + str(wasteSpecies) + "\n"
				
				cmdFileFid.write(strToWrite)
				
			print "\n****\n\n ANALYSIS of the SIMULATION",  str(tmpDir), " IS FINISHED\n"
			
cmdFileFid.close()