#!/usr/bin/python
# -*- coding: latin-1 -*-

import sys, os # Standard librar
import glob
from copy import deepcopy
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
		scc_in_raf = return_scc_in_raf(rctsRAF, cats, rafset[0])
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

def create_chemistry(args, originalSpeciesList, parameters, rctToCat, totCleavage, totCond):
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
		if i % 100 == 0: print "\t\t|- Reaction ", i
		rctType = 1
		if (args.creationMethod == 1) | (args.creationMethod == 4):
			if args.rctRatio > 0: 
				if ran.random() > args.rctRatio: rctType = 0
			else:
				if (totCleavage / (float(totCleavage) + totCond)) <= ran.random(): rctType = 0
		
		# CREATE REACTION 
		# If cleavage or WIM method
		reactionValid = False
		if (rctType == 1) & (nCleavage <= totCleavage):
			rctnew = False
			# Create random cleavage
			while reactionValid == False:
				tmp1, tmp2, tmp3, tmp1id, tmp2id, tmp3id = reactions.createRandomCleavageForCompleteFiringDisk(speciesList, parameters[14], initSpeciesListLength)
				if reactionID > 0:
					if args.revRcts == 0:
						if sum((rcts[:,1] == 1) & (rcts[:,2] == tmp1id) & (rcts[:,3] == tmp3id)) == 0: reactionValid = True
					else:
						reactionValid = True
				else:
					reactionValid = True
				
			# Check if the reaction is new
			if reactionID == 0: rctnew = True
			else: 
				if sum((rcts[:,1] == 1) & (rcts[:,2] == tmp1id) & (rcts[:,3] == tmp2id)) == 0: rctnew = True
				
			# Reaction Structure Creation
			if args.directRctDirection == 1: dirrct = 1
			elif args.directRctDirection == 0: dirrct = 0
			else: dirrct = int(round(ran.random()))
			
			if rctnew: # In the reaction is new
				if reactionID == 0: # If the reaction is the first one
					rcts = np.array([[int(reactionID), int(rctType), tmp1id, tmp2id, tmp3id, int(0), int(239), parameters[34]]])
					reactionID += 1
				else: 
					rcts = np.vstack([rcts,(int(reactionID), int(rctType), tmp1id, tmp2id, tmp3id, int(0), int(239), parameters[34])])	
					reactionID += 1
				nCleavage += 1
					
				if (args.creationMethod == 2) | (args.creationMethod == 4): # If WIM method the reverse reaction is added
					rcts = np.vstack([rcts,(int(reactionID), int(0), tmp1id, tmp2id, tmp3id, int(0), int(239), parameters[33])])	
					reactionID += 1
					nCondensa += 1
			else:
				rct2cat = rcts[(rcts[:,1] == 1) & (rcts[:,2] == tmp1id) & (rcts[:,3] == tmp2id),0]
				if (args.creationMethod == 2) | (args.creationMethod == 4): # If the reverse reaction is present, so the ID is stored
					rct2cat_no_rev = rcts[(rcts[:,1] == 0) & (rcts[:,2] == tmp1id) & (rcts[:,3] == tmp2id),0]
					
					#rct2cat_no_rev = rcts_no_rev[(rcts_no_rev[:,1] == dirrct) & (rcts_no_rev[:,2] == tmp1id) & (rcts_no_rev[:,3] == tmp2id),0]
		else: # condensation
			if (args.creationMethod == 1) | (args.creationMethod == 4):
				rctnew = False
				sub1, sub2, idsub1, idsub2, prod = reactions.createRandomCondensation(speciesList, initSpeciesListLength)
				try:
					tmpprodid = speciesList.index(prod)
				except:
					tmpprodid = len(speciesList)
					speciesList.append(prod)
				
				# The reaction is new if reactinID == 0 (first reaction) or if it is not present 
				if reactionID == 0: rctnew = True
				else:
					if sum((rcts[:,1] == 0) & (rcts[:,3] == idsub1) & (rcts[:,4] == idsub2)) == 0: rctnew = True
				# Reaction Structure Creation
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
					rct2cat = rcts[(rcts[:,1] == 0) & (rcts[:,3] == idsub1) & (rcts[:,4] == idsub2),0]
					if (args.creationMethod == 2) | (args.creationMethod == 4): # If the reverse reaction is present, so the ID is stored
						rct2cat_no_rev = rcts[(rcts[:,1] == 1) & (rcts[:,3] == idsub1) & (rcts[:,4] == idsub2),0]
					
		# A CATALYST IS RANDOMLY (UNIFORM or PREF ATTACHMENT network creation method) ASSIGNED FROM THE SPECIES LIST
		catalyst = -1
		catFound = False
		
		while not catFound:
			if (args.prefAttach == 0) | (i == 0): 
				catalyst = originalSpeciesList.index(ran.choice(originalSpeciesList[len(parameters[14]):initSpeciesListLength-1]))
				weightCat[catalyst] += 1
				pweightCat = [float(i)/sum(weightCat) for i in weightCat]
				
			else:
				#catalyst = choice(range(initSpeciesListLength),p=pweightCat) # TO USE SINCE NUMPY 1.7
				catalyst = species.weightedChoice(pweightCat, range(initSpeciesListLength))
				weightCat[catalyst] += 1
				pweightCat = [float(i)/sum(weightCat) for i in weightCat]
			# IF the selected catalyst is greater than the minimum length and the catalyst does not already catalyze the reaction
			if (len(originalSpeciesList[catalyst]) > args.noCat):
				if rctnew == False:
					if sum((cats[:,1]==catalyst) & (cats[:,2]==rct2cat))==0:
						catFound = True
				else:
					catFound = True
		
		# Forward reaction to catalyze
		if rctnew: 
			if (args.creationMethod == 2) | (args.creationMethod == 4):
				rctsToCat = reactionID - 2
			else:
				rctsToCat = reactionID - 1
		else: 
			rctsToCat = rct2cat
		
		if catalysisID == 0: # if this is the first catalysis
			
			cats = np.array([[int(catalysisID), int(catalyst), int(rctsToCat), int(0), parameters[27], parameters[28], parameters[29], int(1)]])
			catalysisID += 1
				
			if (args.creationMethod == 2) | (args.creationMethod == 4): # IF wim method
				cats = np.vstack([cats,(int(catalysisID), int(catalyst), int(rctsToCat + 1), int(0), parameters[27], parameters[28], parameters[29], int(1))])
				catalysisID += 1

		else: 
			cats = np.vstack([cats,(int(catalysisID), int(catalyst), int(rctsToCat), int(0), parameters[27], parameters[28], parameters[29], int(1))])
			catalysisID += 1
			
			if (args.creationMethod == 2) | (args.creationMethod == 4):
				cats = np.vstack([cats,(int(catalysisID), int(catalyst), int(rctsToCat + 1), int(0), parameters[27], parameters[28], parameters[29], int(1))])
				catalysisID += 1
	
	return rcts, cats, speciesList, rcts_no_rev, cats_no_rev
	
		
	
