#!/usr/bin/python
# -*- coding: latin-1 -*-
'''This python program assess different network structures in term of topological (NON DYNAMICAL) RAF and SCC presence according 
	to different structural parameters. 
'''
import sys, os # Standard librar
import glob
from copy import deepcopy
import numpy as np # Scientific library
import itertools as it
import random as ran
from time import time
from numpy import * 
#from numpy.random import choice  TO USE SINCE NUMPY 1.7
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


#  Input parameters definition 
if __name__ == '__main__':
	parser = ArgumentParser(
				description='This script perform a topological analysis of random catalytic reaction networks.'
				, epilog='''Random Catalytic Reaction Networks Topological Assessment ''') 
	parser.add_argument('-k', '--creationMethod', help='Network creation method (1: Filisetti, 2: Wim, 3: WimNoRevs, 4: Filisetti with revs, DEF: 1)', default='1', type=int)
	parser.add_argument('-R', '--revRcts', help='Reverse reactions are allowed to be created for chance(1: Yes, 0: No, Deafult: No)', default='0', type=int)		
	parser.add_argument('-a', '--prefAttach', help='Type of catalyst choice (1: Preferential Attachment, 0: Random attachment, DEF: 0', default='0', type=int)
	parser.add_argument('-v', '--minavgcon', help='Initial Average Connectivity, level of catalysis, DEF: 0.1', default='0.1', type=float)
	parser.add_argument('-V', '--maxavgcon', help='Final Average Connectivity, level of Catalysis, DEF: 4.0', default='4.0', type=float)
	parser.add_argument('-f', '--lastFood', type=int, help='Last food species ID (default: 5)', default='5')
	parser.add_argument('-o', '--strOut', help='Path for output file storing (Default: ./)', default='./')
	parser.add_argument('-M', '--maxDim', help='Max Dimension of the systems (Default: 10)', default='10', type=int)
	parser.add_argument('-m', '--minDim', help='min Dimension of the systems (Default: 5)', default='5', type=int)
	parser.add_argument('-i', '--iteration', help='How many network instances per dimension are created (Default: 10)', default='10', type=int)
	parser.add_argument('-n', '--noCat', help='Non catalytic max length (default: 2)', default='2', type=int)
	parser.add_argument('-c', '--rctRatio', help='Ratio between cleavages and condensations (default: 0, it means that the actual ratio is used)', default='0', type=float)	
	parser.add_argument('-r', '--randomSeed', help='random seed objects path', default=None)
	parser.add_argument('-d', '--directRctDirection', help='Direction of the forward reaction where necessary (1: cleavage, 0: condensation, 2: random with probability 0.5,Default: 1)', default='1', type=int)

	args = parser.parse_args()
	
	print "\nSystem Initialization..."
		
	# If the architecture of the network allows reverse reaction the revRcts parameter is automatically set to 1
	if (args.creationMethod == 2) | (args.creationMethod == 4): 
		args.revRcts = 1
		print "\n\t!!! WARNING !!! Since the chosen architecture (", args.creationMethod, ") allows the presence of reverse reactions ",\
			" the parameter args.revRcts has been automatically set to 1\n" 
	
	# Create stas folders
	ndn = '_topological_assessment_' +  str(args.creationMethod) + '_' + str(args.lastFood) + '_' \
			+ str(args.noCat) + '_' + str(args.rctRatio) + '_' + str(args.prefAttach)
	newdirAllResults = os.path.join(args.strOut, ndn)
	if not os.path.isdir(newdirAllResults):
		try:
			os.mkdir(newdirAllResults)
		except:
			print "Impossible to create statistic directory", newdirAllResults; sys.exit(1)
	
	# RANDOM SEED INITIALIZATION
	if args.randomSeed == None:
		#! Set randomly the random state and store it in a file
		ran.seed(args.randomSeed)
		writefiles.saveRandomSeed(newdirAllResults)
	else:
		readfiles.loadRandomSeed(args.randomSeed)
	
	print "|- Output Folder: ", args.strOut

	fname_initRafRes = os.path.join(newdirAllResults, '0_initRafAnalysis.csv')
	fname_initRafResSUM = os.path.join(newdirAllResults, '0_initRafAnalysisSUM.csv')
	fname_initRafResLIST = os.path.join(newdirAllResults, '0_initRafAnalysisLIST.csv')
	fname_initRafResALL = os.path.join(newdirAllResults, '0_initRafAnalysisALL.csv')
	fid_initRafRes = open(fname_initRafRes, 'w')
	fid_initRafResSUM = open(fname_initRafResSUM, 'w')
	fid_initRafResLIST = open(fname_initRafResLIST, 'w')
	fid_initRafResALL = open(fname_initRafResALL, 'w')
	strToWrite = "Folder\tP\tAC\tM\tRAFsize\tClosure\tCats\tuRAF\tSCC\tAutoCat\n"
	fid_initRafRes.write(strToWrite)
	strToWrite = "P\tM\tAC\tRAF%\tSCC%\t%SCCinRAF\t%SelfInRAF\tTIME\tTimeCreation\tTimeAnalysis\n"
	fid_initRafResSUM.write(strToWrite)
	
	for maxlength in range(args.minDim,args.maxDim+1): # FOR EACH DIMENSION
		#sys.stdout.flush() # Force save data on file
		increaseYet = True 
		averageConn = args.minavgcon
		while (increaseYet == True) & (averageConn <= args.maxavgcon): # FOR EACH AVERAGE CONNECTIVITY (LEVEL OF CATALYSIS)
			time1 = time()
			raffound = 0
			sccfound = 0
			sccinraffound = 0
			self_sccinraffound = 0
			iterations = int((args.iteration/(maxlength+1-args.minDim)))
			timeCreatinVector = []
			timeAnalysisVector = []
			for instanceID, instance in enumerate((range(iterations))): # FOR EACH ITERATION
				
				# Flush standard output
				sys.stdout.flush()
				
				nCleavage = 0
				nCondensa = 0
				
				# Create the complete population of species according to the Max Length and the alphabet
				alphabet = ['A', 'B']
				speciesList = sp.createCompleteSpeciesPopulation(maxlength, alphabet)
				originalSpeciesList = deepcopy(speciesList)
				
				# compute population carinality
				totSpecies = sp.getTotNumberOfSpeciesFromCompletePop(maxlength)
				
				# Compute overall conceivable number of reactions
				totCleavage = reactions.getNumOfCleavages(speciesList)
				
				if (args.creationMethod == 1) | (args.creationMethod == 4): totCond = reactions.getNumOfCondensations(totSpecies)
				else: totCond = 0
				totRcts = totCleavage + totCond
				
				# If the reaction probability is 0, it is set to the critical value of 1 reaction per species on average
				prob = (1 / float(totRcts)) * averageConn
				
				
				if instanceID == 0: print '|- Max Length: ', maxlength, ' - AVG CON: ', averageConn, ' - Species: ', \
									totSpecies, '- Reactions: ', totRcts, ' <- Condensations: ', totCond, ' - Cleavage: ', totCleavage
				# Compute reactions to catalyse according to the probability
				# rctToCat = int(round(totRcts * totSpecies * prob))
				rctToCat = int(round(totSpecies * averageConn))
				if instanceID == 0: print '\t\t|- P: ', prob, ' , Species (', totSpecies, \
											'), reactions (', totRcts ,'),', rctToCat,\
											' catalysis. Iteration --> ', iterations
				
				initSpeciesListLength = len(speciesList)
				
				conf = (1,1,2000,0,200000,0,0,2,args.lastFood,prob)	
				reactionID = 0
				catalysisID = 0
				#print "\t\t|- NET creation... "
				checkRct = False
				weightCat = [1]*len(speciesList)
				# If the network structure allows the presence of reverse reactions, hence overall reactions must be divided by 2
				if (args.creationMethod == 2) | (args.creationMethod == 4): rctToCat = int(round(rctToCat/2))
				
				pars = [0]*35
				pars[14]="AB";pars[27]=25;pars[28]=50;pars[29]=50
				
				# ARTIFICIAL CHEMISTRY CREATION ----------
				timeNetCreation = time()
				rcts, cats, speciesList, rcts_no_rev, cats_no_rev = network.create_chemistry(args, originalSpeciesList, pars, rctToCat, totCleavage, totCond)
				timeCreatinVector.append(time()-timeNetCreation)
				#print timeCreatinVector
				
				# Create food list
				foodList = range(args.lastFood+1)
				
				# TOPOLOGICAL ANALYSIS
				#print "\t\t|- RAF searching step..."
				timeAnalysis = time()
				netouts = network.net_analysis_of_static_graphs(fid_initRafRes, fid_initRafResALL, fid_initRafResLIST, 'tmpDir', prob, averageConn, rcts, cats, foodList, maxlength)
				timeAnalysisVector.append(time() - timeAnalysis)
				#print timeAnalysisVector

				# If RAFs are present, the presence of possible SCCs in RAF is assessed
				if len(netouts[0][2]) > 0: 
					raffound += 1
					rctsRAF = rcts[np.any(rcts[:, 0] == np.expand_dims(netouts[0][2],1), 0), :]
					scc_in_raf = network.return_scc_in_raf(rctsRAF, netouts[0][0], cats)
					if scc_in_raf[1] > 0:
						sccinraffound += 1
					if scc_in_raf[2] > 0:
						self_sccinraffound += 1
				if netouts[1][4] > 0: sccfound += 1

				del rcts
				del cats
		
			time2 = time()
			timeAvgCreation = sum(timeCreatinVector) / float(iterations)
			timeAvgAnalysis = sum(timeAnalysisVector) / float(iterations)
			percRAFFounds = raffound/float(iterations)
			percSccFounds = sccfound/float(iterations)
			if raffound > 0: percScc_in_rafFounds = sccinraffound / float(raffound)
			else: percScc_in_rafFounds = 0
			if raffound > 0: perc_self_Scc_in_rafFounds = self_sccinraffound / float(raffound)
			else: perc_self_Scc_in_rafFounds = 0
			#print percScc_in_rafFounds
			#print perc_self_Scc_in_rafFounds
			#raw_input("ciao")
			
			timet = time2-time1
			print "\t\t\t %RAF *** ", percRAFFounds, " *** %SCC ***", percSccFounds," *** TIME: ", time2 - time1, " Cleavages: ", nCleavage, " - Condensations: ", nCondensa
			if percRAFFounds >= 0.99: 
				print "\t\t\t Max number of RAFs found"
				increaseYet = False # If RAF is always found so next average conn is assessed
			strToWrite = str(prob) + "\t" + str(maxlength) + "\t" + str(averageConn) + "\t" + str(percRAFFounds) + "\t" + str(percSccFounds) + \
						"\t" + str(percScc_in_rafFounds) + "\t" + str(perc_self_Scc_in_rafFounds) + "\t" + str(timet) + "\t" + str(timeAvgCreation) + \
						"\t" + str(timeAvgAnalysis) + "\n"
			fid_initRafResSUM.write(strToWrite)
			averageConn += 0.1
	# Close text files
	fid_initRafRes.close()
	fid_initRafResLIST.close()
	fid_initRafResALL.close()
	fid_initRafResSUM.close()	
	
		