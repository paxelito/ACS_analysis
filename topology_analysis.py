#!/usr/bin/python
# -*- coding: latin-1 -*-
'''Script to order the analysis of the divergences in time. 
'''
import sys, os # Standard librar
import glob
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


#  Input parameters definition 
if __name__ == '__main__':
	parser = ArgumentParser(
				description='This script perform a topological analysis of random catalytic reaction networks.'
				, epilog='''Random Catalytic Reaction Networks Topological Assessment ''') 
	parser.add_argument('-k', '--creationMethod', help='Network creation method (1: Filisetti, 2: Wim, 3: WimNoRevs, 4: Filisetti with revs, DEF: 1)', default='1', type=int)
	parser.add_argument('-f', '--lastFood', type=int, help='Last food species ID (deafult: 5)', default='5')
	parser.add_argument('-o', '--strOut', help='Path for output file storing (Default: ./)', default='./')
	parser.add_argument('-M', '--maxDim', help='Max Dimension of the systems (Default: 10)', default='10', type=int)
	parser.add_argument('-m', '--minDim', help='min Dimension of the systems (Default: 5)', default='5', type=int)
	parser.add_argument('-i', '--iteration', help='How many network instances per dimension are created (Default: 10)', default='10', type=int)
	parser.add_argument('-n', '--noCat', help='Non catalytic max length (default: 2)', default='2', type=int)
	parser.add_argument('-c', '--rctRatio', help='Ratio between cleavages and condensations (default: 0, it means that the actual ratio is used)', default='0', type=float)	
	parser.add_argument('-r', '--randomSeed', help='random seed', type=int, default=None)
	args = parser.parse_args()
	
	print "\nSystem Initialization..."
	ran.seed(args.randomSeed)
	
	# Create stas folders
	ndn = '_0_new_allStatResults_' +  str(args.creationMethod) + '_' + str(args.lastFood) + '_' + str(args.noCat)
	newdirAllResults = os.path.join(args.strOut, ndn)
	if not os.path.isdir(newdirAllResults):
		try:
			os.mkdir(newdirAllResults)
		except:
			print "Impossible to create statistic directory", newdirAllResults; sys.exit(1)
	
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
	strToWrite = "P\tM\tAC\tRAF%\tSCC%\t%SCCinRAF\t%SelfInRAF\tTIME\n"
	fid_initRafResSUM.write(strToWrite)
	
	for maxlength in range(args.minDim,args.maxDim+1): 
		#sys.stdout.flush() # Force save data on file
		increaseYet = True 
		averageConn = 0.1
		while (increaseYet == True) & (averageConn < 3):
			time1 = time()
			raffound = 0
			sccfound = 0
			sccinraffound = 0
			self_sccinraffound = 0
			iterations = int((args.iteration/(maxlength+1-args.minDim)))
			for instanceID, instance in enumerate((range(iterations))): 
				nCleavage = 0
				nCondensa = 0
				
				# Create the complete population of species according to the Max Length and the alphabet
				alphabet = ['A', 'B']
				species = sp.createCompleteSpeciesPopulation(maxlength, alphabet)
				
				# compute population carinality
				totSpecies = sp.getTotNumberOfSpeciesFromCompletePop(maxlength)
				
				# Compute overall conceivable number of reactions
				totCleavage = reactions.getNumOfCleavages(species)
				
				if args.creationMethod == 1: totCond = reactions.getNumOfCondensations(totSpecies)
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
											'), , reactions (', totRcts ,'),', rctToCat,\
											' catalysis. Iteration --> ', iterations
				
				initSpeciesListLength = len(species)
				
				conf = (1,1,2000,0,200000,0,0,2,args.lastFood,prob)	
				reactionID = 0
				catalysisID = 0
				#print "\t\t|- NET creation... "
				checkRct = False
				for i in range(rctToCat):
					# Select if condensation of cleavage according to the total number of reactions
					#if rctToCat > 1000: 
					#	if i % 1000 == 0: print "\t\t\t|- Reaction ", i, " - species list L: ", len(species)
					rctType = 1
					if (args.creationMethod == 1) | (args.creationMethod == 4):
						if args.rctRatio > 0: 
							if ran.random() > args.rctRatio: rctType = 0
						else:
							if (totCleavage / (float(totCleavage) + totCond)) <= ran.random(): rctType = 0
					
					 
					# If cleavage
					if (rctType == 1) & (nCleavage <= totCleavage):
						# Select species to cleave
						rctnew = False
						# Create random cleavage
						molToCleav, tmp1, tmp2, tmp1id, tmp2id, find1 = reactions.createRandomCleavage(species, alphabet, initSpeciesListLength)
						# Check if the reaction is new
						if i == 0: rctnew = True
						else:
							if find1 == True:
								if sum((rcts[:,1] == 1) & (rcts[:,2] == species.index(molToCleav)) \
									& (rcts[:,3] == tmp1id)) == 0: rctnew = True
							
						
						if rctnew:
							if reactionID == 0:
								rcts = np.array([[int(reactionID), int(rctType), species.index(molToCleav), tmp1id, tmp2id, int(0), int(0), int(0)]])
								reactionID += 1
								nCleavage += 1
							else: 
								rcts = np.vstack([rcts,(int(reactionID), int(rctType), species.index(molToCleav), tmp1id, tmp2id, int(0), int(0), int(0))])	
								reactionID += 1
								nCleavage += 1
								
							if (args.creationMethod == 2) | (args.creationMethod == 4): # if reverse reaction are allowed (methods 2 and 4 :: WIM and OUR)
								rcts = np.vstack([rcts,(int(reactionID), int(0), species.index(molToCleav), tmp1id, tmp2id, int(0), int(0), int(0))])	
								reactionID += 1
								nCondensa += 1
						else:
							rct2cat = rcts[(rcts[:,1] == 1) & (rcts[:,2] == species.index(molToCleav)) & (rcts[:,3] == tmp1id),0]
					else: # condensation
						if (args.creationMethod == 1) | (args.creationMethod == 4):
							rctnew = False
							sub1, sub2, idsub1, idsub2, prod = reactions.createRandomCondensation(species, initSpeciesListLength)
							try:
								tmpprodid = species.index(prod)
							except:
								tmpprodid = len(species)
								species.append(prod)
							
							if i == 0: rctnew = True
							else:
								if sum((rcts[:,1] == 0) & (rcts[:,3] == idsub1) & (rcts[:,4] == idsub2)) == 0: rctnew = True
							
							# Reaction Structure Creation
							if rctnew:
								if reactionID == 0: 
									rcts = np.array([[int(reactionID), int(rctType), tmpprodid, idsub1, idsub2, int(0), int(0), int(0)]])
									reactionID += 1
								else: 
									rcts = np.vstack([rcts,(int(reactionID), int(rctType), tmpprodid, idsub1, idsub2, int(0), int(0), int(0))])	
									reactionID += 1
								nCondensa += 1
							else:
								rct2cat = rcts[(rcts[:,1] == 0) & (rcts[:,3] == idsub1) & (rcts[:,4] == idsub2),0]
								
							if args.creationMethod == 4: # if reverse reaction are allowed (methods 2 and 4 :: WIM and OUR)
								rcts = np.vstack([rcts,(int(reactionID), int(1), tmpprodid, idsub1, idsub2, int(0), int(0), int(0))])	
								reactionID += 1
								nCondensa += 1
					# A CATALYST IS RANDOMLY ASSIGNED FROM THE SPECIES LIST
					catalyst = -1
					catFound = False
					
					while not catFound: 
						catalyst = species.index(ran.choice(species[len(alphabet):]))
						if (len(species[catalyst]) > args.noCat):
							if rctnew == False:
								if sum((cats[:,1]==catalyst) & (cats[:,2]==rct2cat))==0:
									catFound = True
							else:
								catFound = True
						
					if rctnew: rctsToCat = reactionID - 1
					else: rctsToCat = rct2cat
					
					if (args.creationMethod == 2) | (args.creationMethod == 4): rctsToCat = reactionID - 2
					if catalysisID == 0: 
						cats = np.array([[int(catalysisID), int(catalyst), int(rctsToCat), int(0), float(0.5), float(0.25), float(0.5), ran.randint(1,2)]])
						catalysisID += 1
						if (args.creationMethod == 2) | (args.creationMethod == 4): # IF wim method
							cats = np.vstack([cats,(int(catalysisID), int(catalyst), int(rctsToCat + 1), int(0), float(0.5), float(0.25), float(0.5), ran.randint(1,2))])
							catalysisID += 1
					else: 
						cats = np.vstack([cats,(int(catalysisID), int(catalyst), int(rctsToCat), int(0), float(0.5), float(0.25), float(0.5), ran.randint(1,2))])
						catalysisID += 1
						if (args.creationMethod == 2) | (args.creationMethod == 4):
							cats = np.vstack([cats,(int(catalysisID), int(catalyst), int(rctsToCat + 1), int(0), float(0.5), float(0.25), float(0.5), ran.randint(1,2))])
							catalysisID += 1
				#print rcts
				#print cats
				#print cats.shape
				#print rcts.shape
				#raw_input("ca0")
				# Create food list
				foodList = range(args.lastFood+1)
				#print "\t\t|- RAF searching step..."
				netouts = network.net_analysis_of_static_graphs(fid_initRafRes, fid_initRafResALL, fid_initRafResLIST, 'tmpDir', prob, averageConn, rcts, cats, foodList, maxlength)
				#print netouts
				if len(netouts[0][2]) > 0: 
					raffound += 1
					rctsRAF = rcts[np.any(rcts[:, 0] == np.expand_dims(netouts[0][2],1), 0), :]
					scc_in_raf = network.return_scc_in_raf(rctsRAF, cats, netouts[0][0])
					if scc_in_raf[1] > 0:
						sccinraffound += 1
					if scc_in_raf[2] > 0:
						self_sccinraffound += 1
				if netouts[1][4] > 0: sccfound += 1
				#print raffound
				#raw_input("ciao")

				del rcts
				del cats
		
			time2 = time()
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
						"\t" + str(percScc_in_rafFounds) + "\t" + str(perc_self_Scc_in_rafFounds) + "\t" + str(timet) + "\n"
			fid_initRafResSUM.write(strToWrite)
			averageConn += 0.1
	# Close text files
	fid_initRafRes.close()
	fid_initRafResLIST.close()
	fid_initRafResALL.close()
	fid_initRafResSUM.close()	
	
		