#!/usr/bin/python
# -*- coding: latin-1 -*-

import sys, os # Standard librar
import glob
from copy import deepcopy
from time import time
import random as ran
import numpy as np # Scientific library
from numpy import * 
#from numpy.random import choice  TO USE SINCE NUMPY 1.7
from argparse import ArgumentParser
try:
    from pylab import *
except:
    pass
   
from ..graph import raf
from ..graph import scc
from ..dyn import dynamics as dn
from ..IO import writefiles
from ..model import reactions
from lib.model import species
   
   
def removeRareRcts(graph, dt, life, nrg, deltat):
	if shape(graph)[0] > 1:
		graph[:,life] = graph[:,life] + deltat
		graph[:,nrg] = graph[:,dt] - graph[:,life]
		graph = graph[graph[:,nrg]>0,:]
		return graph
	else:
		return graph

def fixCondensationReaction(m1, m2, m3, rcts):
	
	if sum((rcts[:,2] == m1) & (rcts[:,3] == m2) & (rcts[:,4] == m3)) >= 1:
		return m1, m2, m3
	elif sum((rcts[:,2] == m1) & (rcts[:,3] == m3) & (rcts[:,4] == m2)) >= 1:
		return m1, m3, m2
	else: 
		print m1, m2, m3
		print "ERROR!!!!"
		sys.exit(1)

# BRIDGE FUNCTION TO DETECT RAFs in INITIAL SETS
def net_analysis_of_static_graphs(fid_initRafRes, fid_initRafResALL, fid_initRafResLIST, tmpDir, rctProb, avgCon, rcts, cats, foodList, maxDim,debug=False):
	rafset = raf.rafsearch(rcts, cats, foodList,debug) #RAF search 
	stdgraph = scc.createNetXGraph(rcts,cats)
	sccsets = scc.diGraph_netX_stats(stdgraph)
	ErctP = "%.4g" % rctProb
	strToWrite = tmpDir + "\t" + str(ErctP) + "\t" + str(avgCon) + "\t" + str(maxDim) + "\t" + \
				 str(len(rafset[2])) + "\t" + str(len(rafset[0])) + "\t" + str(len(rafset[3])) + \
				 "\t" + str(rafset[4]) + "\t" + str(sccsets[4]) + "\t" + str(sccsets[2]) + "\n"
	fid_initRafRes.write(strToWrite)
	writefiles.write_init_raf_list(fid_initRafResLIST, rafset, tmpDir)
	writefiles.write_init_raf_all(fid_initRafResALL, rafset, tmpDir, rcts, cats)
	return rafset, sccsets

# BRIDGE FUNCTION TO DETECT RAFs in DYNAMICS
def net_analysis_of_dynamic_graphs(fid_dynRafRes, tmpTime, rcts, cats, foodList, growth=False, rctsALL=None, catsALL=None, completeRCTS=None,debug=False):
	#print rcts
	#print cats
	rafset = raf.rafsearch(rcts, cats, foodList,debug) #RAF search
	stdgraph = scc.createNetXGraph(rcts,cats) #netX graph creation
	sccsets = scc.diGraph_netX_stats(stdgraph) # SCC analysis
	if growth == True: 
		rafsetALL = raf.rafsearch(rctsALL, catsALL, foodList) #RAF search
		stdgraphALL = scc.createNetXGraph(rctsALL,catsALL) # netX graph creation
		sccsetsALL = scc.diGraph_netX_stats(stdgraphALL) # SCC analysis
		
	strRAF = '' 
	# If RAF analysis is made in dynamical temporary structures a trnaslation in real net must be done
	if len(rafset[2]) > 0: 
		rctsRAF = rcts[np.any(rcts[:, 0] == np.expand_dims(rafset[2],1), 0), :]
		scc_in_raf = return_scc_in_raf(rctsRAF, rafset[0], cats)
	else:
		scc_in_raf = 0,0,0,0,0,0
	 
	if completeRCTS != None: convRAF = raf.findRAFrcts(rafset[2],rcts,completeRCTS)
	else: convRAF = rafset[2]
	if len(convRAF) > 0: 		
		for x in convRAF: strRAF += str(x) + '\t'	
	
	#t	#CL	#RAF	#SCC	#SelfCats	RAF		
		
	if growth == False: strToWrite = str(tmpTime) + "\t" + str(len(rafset[0])) + "\t" + str(rafset[4]) + "\t" + str(sccsets[4]) + "\t" + str(sccsets[2]) + "\t" + str(scc_in_raf[1])+ "\t" + strRAF + "\n"
	else: strToWrite = str(tmpTime) + "\t" + str(len(rafset[0])) + "\t" + str(rafset[4]) + str(len(rafsetALL[0])) + "\t" + str(rafsetALL[4]) + "\t" + str(sccsets[4]) + "\t" + str(sccsets[2]) + "\t" + str(scc_in_raf[1])+ "\t" + strRAF + "\n"
	fid_dynRafRes.write(strToWrite)
	return rafset

def return_scc_in_raf(tmpRAF, tmpClosure, tmpCats):
	stdgraph = scc.createNetXGraphForRAF(tmpRAF, tmpClosure, tmpCats)
	sccsets = scc.diGraph_netX_stats(stdgraph)
	return sccsets 

def create_chemistry(args, originalSpeciesList, parameters, rctToCat, totCleavage, totCond, tmpac, autocat=True):
	#print '\t|- Create chemistry'
	speciesList = deepcopy(originalSpeciesList)
	initSpeciesListLength = len(speciesList) # Initial cardinality of the species list (to avoid recursive multiple species evaluation)
	reactionID = 0
	catalysisID = 0
	catalysisID_no_rev = 0
	reactionID_no_rev = 0
	nCleavage = 0
	nCondensa = 0
	rcts_no_rev = []
	cats_no_rev = []
	weightCat = [1]*len(originalSpeciesList)
	#print "\t\t|- NET creation... "
	checkRct = False
	
	for i in range(rctToCat):
		if i > 0:
			if i % 1000 == 0: print '\t|- Reaction ', i+1
		rctTime = time()
		rctType = 1
		if (args.creationMethod == 1) | (args.creationMethod == 4):
			if args.rctRatio < 1: 
				if ran.random() > args.rctRatio: rctType = 0
			else:
				if (totCleavage / (float(totCleavage) + totCond)) <= ran.random(): rctType = 0
		
		# CREATE REACTION 
		# If cleavage or WIM method
		reactionValid = False
		if (rctType == 1) & (nCleavage <= totCleavage):
			rctnew = False
			revType = 0 # By default the reverse reaction of the cleavage is condensation
			# Create random cleavage (or 50% cleavage/condensation in case of creation method 3 (WIM without reverse reactions)
			# If the creation method is 3 (WIM without reverse reactions) the nature of the reaction is randomly selected
			# according to the rctRatio parameter (0 means all cleavages) 
			if (args.creationMethod == 3):
				if args.rctRatio < 1:
					if ran.random() < args.rctRatio:
						rctType = 1
						revType = 0
					else:
						rctType = 0
						revType = 1				
			
			# Reaction Structure Creation	
			while reactionValid == False:
				tmp1, tmp2, tmp3, tmp1id, tmp2id, tmp3id = reactions.createRandomCleavageForCompleteFiringDisk(speciesList, parameters[14], initSpeciesListLength)
				if reactionID > 0:
					if args.revRcts == 0: # If rev rcts are neglected, reverse reactions already present are searched
						if sum((rcts[:,1] == revType) & (rcts[:,2] == tmp1id) & (rcts[:,3] == tmp2id)) == 0: reactionValid = True
					else:
						reactionValid = True
				else:
					reactionValid = True
				
			# Check if the reaction is new
			if reactionID == 0: rctnew = True
			else: 
				if sum((rcts[:,1] == rctType) & (rcts[:,2] == tmp1id) & (rcts[:,3] == tmp2id)) == 0: rctnew = True
			
			if rctnew: # In the reaction is new
				if reactionID == 0: # If the reaction is the first one
					rcts = np.array([[int(reactionID), int(rctType), tmp1id, tmp2id, tmp3id, int(0), int(239), parameters[34]]])
					reactionID += 1
				else: 
					rcts = np.vstack([rcts,(int(reactionID), int(rctType), tmp1id, tmp2id, tmp3id, int(0), int(239), parameters[34])])	
					reactionID += 1
				if rctType == 1: nCleavage += 1
				else: nCondensa += 1
					
				if (args.creationMethod == 2) | (args.creationMethod == 4): # If WIM method the reverse reaction is added
					rcts = np.vstack([rcts,(int(reactionID), int(revType), tmp1id, tmp2id, tmp3id, int(0), int(239), parameters[33])])	
					reactionID += 1
					if revType == 1: nCleavage += 1
					else: nCondensa += 1
			else:
				rct2cat = rcts[(rcts[:,1] == rctType) & (rcts[:,2] == tmp1id) & (rcts[:,3] == tmp2id),0]
				if (args.creationMethod == 2) | (args.creationMethod == 4): # If the reverse reaction is present, so the ID is stored
					if rct2cat % 2 == 1: rct2cat -= 1 # If the selected reaction is already present and the ID is odd it means that the "direct" reaction is the previous one. 
					#rct2cat_no_rev = rcts[(rcts[:,1] == 0) & (rcts[:,2] == tmp1id) & (rcts[:,3] == tmp2id),0]	
		
		else: # condensation
			if (args.creationMethod == 1) | (args.creationMethod == 4):
				rctnew = False
				while reactionValid == False:
					sub1, sub2, idsub1, idsub2, prod = reactions.createRandomCondensation(speciesList, initSpeciesListLength)
					if reactionID > 0:
						if args.revRcts == 0: # If rev rcts are neglected, reverse reactions already present are searched
							if sum((rcts[:,1] == 1) & (rcts[:,3] == idsub1) & (rcts[:,4] == idsub2)) == 0: reactionValid = True
						else:
							reactionValid = True
					else:
						reactionValid = True
				
				# If the reaction is valid, the ID product is selected or created according with the novelty of the product 	
				try:
					tmpprodid = speciesList.index(prod)
				except:
					tmpprodid = len(speciesList)
					speciesList.append(prod)
				
				# The reaction is new if reactinID == 0 (first reaction) or if it is not present 
				if reactionID == 0: rctnew = True
				else:
					if sum((rcts[:,1] == 0) & (rcts[:,3] == idsub1) & (rcts[:,4] == idsub2)) == 0: rctnew = True
					
				# New Reaction is now added (np.vstack) to the creation list (if the list is empty, so the list will be created (np.array))
				if rctnew:
					if reactionID == 0: 
						rcts = np.array([[int(reactionID), int(rctType), tmpprodid, idsub1, idsub2, int(0), int(239), parameters[33]]])
						reactionID += 1
					else: 
						rcts = np.vstack([rcts,(int(reactionID), int(rctType), tmpprodid, idsub1, idsub2, int(0), int(239), parameters[33])])	
						reactionID += 1
					nCondensa += 1
					
					 # if reverse reaction are present (methods 2 and 4 :: WIM and FILISETTI with reverse reactions)
					if args.creationMethod == 4:
						rcts = np.vstack([rcts,(int(reactionID), int(1), tmpprodid, idsub1, idsub2, int(0), int(239), parameters[34])])	
						reactionID += 1
						nCleavage += 1
				else:
					# IF the reaction is now new, the ID is selected in order to set the catalysis
					rct2cat = rcts[(rcts[:,1] == 0) & (rcts[:,3] == idsub1) & (rcts[:,4] == idsub2),0]
					if (args.creationMethod == 2) | (args.creationMethod == 4): # If the reverse reaction is present, so the ID is stored
						if rct2cat % 2 == 1: rct2cat -= 1 # If the selected reaction is already present and the ID is odd it means that the "direct" reaction is the previous one. 
						#rct2cat_no_rev = rcts[(rcts[:,1] == 1) & (rcts[:,3] == idsub1) & (rcts[:,4] == idsub2),0]
					
		# A CATALYST IS RANDOMLY (UNIFORM or PREF ATTACHMENT network creation method) ASSIGNED FROM THE SPECIES LIST
		catalyst = -1
		catFound = False # Flag variable checking the validity of the catalyst
		noSameRct = False # Flag variable to check if the catalyst has already catalysed a given reaction
		
		#if round(tmpac,1) == 1.1: print rcts
				
		while not noSameRct:
			
			# Define reaction ID
			if rctnew: 
				if (args.creationMethod == 2) | (args.creationMethod == 4):
					rctsToCat = reactionID - 2
				else:
					rctsToCat = reactionID - 1
				
				noSameRct = True # correct catalysis, since the reaction is new, it is impossible that it has been already catalysed
				
			else: 
				rctsToCat = rct2cat
				
			# Search for the catalyst 
			while not catFound:
				if (args.prefAttach == 0) | (i == 0): 
					catalyst = originalSpeciesList.index(ran.choice(originalSpeciesList[len(parameters[14]):initSpeciesListLength-1]))
					if (args.prefAttach == 1):
						weightCat[catalyst] += 1.00
						#pweightCat = [float(i)/sum(weightCat) for i in weightCat]
						#pweightCat = map(lambda x: float(x)/sum(weightCat), weightCat)
						np_weightCat = np.array(weightCat)
						pweightCat = np_weightCat / sum(np_weightCat)
					
				else:
					#catalyst = choice(range(initSpeciesListLength),p=pweightCat) # TO USE SINCE NUMPY 1.7
					catalyst = species.weightedChoice(pweightCat, range(initSpeciesListLength))
					weightCat[catalyst] += 1.00
					
					#pweightCat = [float(i)/sum(weightCat) for i in weightCat]
					#pweightCat = map(lambda x: float(x)/sum(weightCat), weightCat)
					np_weightCat = np.array(weightCat)
					pweightCat = np_weightCat / sum(np_weightCat)
						
				# IF the selected catalyst is greater than the minimum length and the catalyst does not already catalyze the reaction
				if (len(originalSpeciesList[catalyst]) > args.noCat):
					if rctnew == False:
						if sum((cats[:,1]==catalyst) & (cats[:,2]==rct2cat))==0:
							if not autocat:
								if int(rctType) == 1:
									if (catalyst == int(rcts[rctsToCat,3])) or (catalyst == int(rcts[rctsToCat,4])): catFound = False
									else: catFound = True
								else:
									if (catalyst == int(rcts[rctsToCat,2])): catFound = False
									else: catFound = True								
							else:
								catFound = True
					else:
						if not autocat:
							if int(rctType) == 1:
								if (catalyst == int(rcts[rctsToCat,3])) or (catalyst == int(rcts[rctsToCat,4])): catFound = False
								else: catFound = True
							else:
								if (catalyst == int(rcts[rctsToCat,2])): catFound = False
								else: catFound = True								
						else:
							catFound = True
						
			if not autocat:
				if int(rctType) == 1:
					if (catalyst == int(rcts[rctsToCat,3])) or (catalyst == int(rcts[rctsToCat,4])): 
						print catalyst, " || ", int(rcts[rctsToCat,2]), " -> ", int(rcts[rctsToCat,3]), " + ", int(rcts[rctsToCat,4])
						raw_input("W")
				else:
					if (catalyst == int(rcts[rctsToCat,2])): 
						print catalyst, " || ", int(rcts[rctsToCat,2]), " <- ", int(rcts[rctsToCat,3]), " + ", int(rcts[rctsToCat,4])
						raw_input("W")
			
			# Store catalysis
			if catalysisID == 0: # if this is the first catalysis
				
				cats = np.array([[int(catalysisID), int(catalyst), int(rctsToCat), int(0), parameters[27], parameters[28], parameters[29], int(1)]])
				catalysisID += 1
					
				if (args.creationMethod == 2) | (args.creationMethod == 4): # IF wim method
					cats = np.vstack([cats,(int(catalysisID), int(catalyst), int(rctsToCat + 1), int(0), parameters[27], parameters[28], parameters[29], int(1))])
					catalysisID += 1
				
				noSameRct = True # correct catalysis, since the catalysis is new, it is impossible that it has been already catalysed	
	
			else: 
				
				# If the double catalysis must be assessed
				if not noSameRct: 
					if sum((cats[:,1]==int(int(catalyst)))&(cats[:,2]==int(rctsToCat))) == 0: # IF the reaction has not yet been catalysed by the catalyst
						cats = np.vstack([cats,(int(catalysisID), int(catalyst), int(rctsToCat), int(0), parameters[27], parameters[28], parameters[29], int(1))])
						catalysisID += 1
						if (args.creationMethod == 2) | (args.creationMethod == 4):
							cats = np.vstack([cats,(int(catalysisID), int(catalyst), int(rctsToCat + 1), int(0), parameters[27], parameters[28], parameters[29], int(1))])
							catalysisID += 1
						noSameRct = True # Reaction definitely valid!!!
					else:
						catFound = False # New catalyst must be found
						print "\t\t\t WARNING: The catalyst ", int(catalyst), " already catalyses the reaction ", int(rctsToCat), ", new catalyst must be found"	
				else: # if the catalysis is valid
					cats = np.vstack([cats,(int(catalysisID), int(catalyst), int(rctsToCat), int(0), parameters[27], parameters[28], parameters[29], int(1))])
					catalysisID += 1
					if (args.creationMethod == 2) | (args.creationMethod == 4):
						cats = np.vstack([cats,(int(catalysisID), int(catalyst), int(rctsToCat + 1), int(0), parameters[27], parameters[28], parameters[29], int(1))])
						catalysisID += 1
		
		#if round(tmpac,1) == 1.1: 
		#	print cats
		#	raw_input("ss")				
	
	return rcts, cats, speciesList, rcts_no_rev, cats_no_rev
	
		
	
