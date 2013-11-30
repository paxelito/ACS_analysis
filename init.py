	#!/usr/bin/env python
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
from lib.dyn import dynamics as dn


#ÊInput parameters definition 
if __name__ == '__main__':
	parser = ArgumentParser(
				description='This script re-arrange results in a more friendly way from the angle analysis in time.'
				, epilog='''File with angle trajectories are created. ''') 
	parser.add_argument('-a', '--sysType', help='System Architecture [1:CLOSE, 2:PROTO, 3:CSTR], deafult: 1', default='1')
	parser.add_argument('-k', '--creationMethod', help='Network creation method (1: Filisetti, 2: Wim, DEF: 1)', default='1', type=int)
	parser.add_argument('-f', '--lastFood', type=int, help='Last food species ID (deafult: 5)', default='5')
	parser.add_argument('-s', '--strPath', help='Path where files are stored', default='./')
	parser.add_argument('-o', '--strOut', help='Path for output file storing', default='./')
	parser.add_argument('-M', '--maxDim', help='Max Dimension of the systems (Default: 10)', default='10', type=int)
	parser.add_argument('-m', '--minDim', help='min Dimension of the systems (Default: 5)', default='5', type=int)
	parser.add_argument('-i', '--iteration', help='How many network instances per dimension are created (Default: 10)', default='10', type=int)
	parser.add_argument('-p', '--rctProb', help='Reaction Probability, if 0 args.avgCon will be used', default='0', type=float)
	parser.add_argument('-v', '--avgCon', help='Tuple containing the average connectivity range (deafult: [0.5, 0.05, 1])	', nargs=3, metavar=('fr', 'st', 'en'), default= [0.5, 0.05, 1])
	parser.add_argument('-n', '--noCat', help='Non catalytic max length (default: 2)', default='2', type=float)
	parser.add_argument('-c', '--rctRatio', help='Ratio between cleavages and condensations (default: 0, it means that the actual ratio is used)', default='0', type=float)	
	parser.add_argument('-r', '--randomSeed', help='random seed', type=int, default=None)
	args = parser.parse_args()
	
	print "\nSystem Initialization..."
	ran.seed(args.randomSeed)
	
	# Create stas folders
	ndn = '_0_new_allStatResults'
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
	strToWrite = "Folder\tP\tAC\tM\tRAFsize\tClosure\tCats\tuRAF\n"
	fid_initRafRes.write(strToWrite)
	strToWrite = "P\tM\tAC\tRAF%\tTIME\n"
	fid_initRafResSUM.write(strToWrite)
	
	for maxlength in range(args.minDim,args.maxDim+1): # For each dimension
		#avgCon = dn.rangeFloat(float(args.avgCon[0]), float(args.avgCon[1]), float(args.avgCon[2]))
		#sys.stdout.flush() # Force save data on file
		if args.rctProb == 0: # if reaction probability is 0
			increaseYet = True 
			averageConn = 0.1
			while increaseYet == True:
				time1 = time()
				raffound = 0
				iterations = int((args.iteration/(maxlength+1-args.minDim)))
				for instanceID, instance in enumerate((range(iterations))): # For each instance of network
					nCleavage = 0
					nCondensa = 0
					# Create all the the species starting from the alphabet
					alphabet = ['A', 'B']
					species = []
					for i in range(maxlength): species.extend(map(''.join,it.product(alphabet, repeat=i+1)))
					
					# compute number of cleavage
					totSpecies = 2 ** (maxlength + 1) - 2
					# Compute overall conceivable number of reactions
					totCleavage = sum(map(lambda x: len(x)-1,species))
					if args.creationMethod == 1:
						totCond = totSpecies ** 2
						totRcts = totCleavage + totCond
					else:
						totRcts = totCleavage
						totCond = 0
					
					# If the reaction probability is 0, it is set to the critical value of 1 reaction per species on average
					if args.rctProb == 0: prob = (1 / float(totRcts)) * averageConn
					else: prob = args.rctProb
					
					if instanceID == 0: print '|- Max Length: ', maxlength, ' - AVG CON: ', averageConn, ' - Species: ', \
										totSpecies, '- Reactions: ', totRcts, ' <- Condensations: ', totCond, ' - Cleavage: ', totCleavage
					# Compute reactions to catalyse according to the probability
					rctToCat = int(round(totRcts * totSpecies * prob))
					if instanceID == 0: print '\t\t|- According to probability ', prob, ' and the number of species (', totSpecies,'), ', \
										rctToCat, ' catalysis will be created. Itertion --> ', iterations
					
					initSpeciesListLength = len(species)
					
					#return (ngen,nsim,totTimes,nrgType,totalRcts,nrgConc,influx_rate,maxLOut,lfdsID,rctProb)
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
						if args.creationMethod == 1:
							if args.rctRatio > 0: 
								if ran.random() > args.rctRatio: rctType = 0
							else:
								if (totCleavage / (float(totCleavage) + totCond)) <= ran.random(): rctType = 0
						
						# CREATE REACTION 
						# If cleavage
						if (rctType == 1) & (nCleavage <= totCleavage):
							# Select species to cleave
							rctnew = False
							molToCleav = ran.choice(species[len(alphabet):initSpeciesListLength-1])
							cutPt = ran.randint(1,len(molToCleav)-1) 
							tmp1 = molToCleav[0:cutPt]
						
							try: 
								tmp1id = species.index(tmp1)
								find1 = True
							except: 
								tmp1id = len(species)+1
								find1 = False
								
							tmp2 = molToCleav[cutPt:len(molToCleav)]
							try: 
								tmp2id = species.index(tmp2)
							except: 
								if find1 == True: tmp2id = len(species) + 1   
								else: tmp2id = len(species)+2
							
							if i == 0: rctnew = True
							else:
								if find1 == True:
									if sum((rcts[:,1] == 1) & (rcts[:,2] == species.index(molToCleav)) \
										& (rcts[:,3] == tmp1id)) == 0: rctnew = True
								
							# Reaction Structure Creation
							if rctnew:
								if reactionID == 0:
									rcts = np.array([[int(reactionID), int(rctType), species.index(molToCleav), tmp1id, tmp2id, int(0), int(0), int(0)]])
									reactionID += 1
									nCleavage += 1
								else: 
									rcts = np.vstack([rcts,(int(reactionID), int(rctType), species.index(molToCleav), tmp1id, tmp2id, int(0), int(0), int(0))])	
									reactionID += 1
									nCleavage += 1
									
								if args.creationMethod == 2: # If WIM method the reverse reaction is added
									rcts = np.vstack([rcts,(int(reactionID), int(0), species.index(molToCleav), tmp1id, tmp2id, int(0), int(0), int(0))])	
									reactionID += 1
									nCondensa += 1
							else:
								rct2cat = rcts[(rcts[:,1] == 1) & (rcts[:,2] == species.index(molToCleav)) & (rcts[:,3] == tmp1id),0]
						else: # condensation
							if args.creationMethod == 1:
								rctnew = False
								sub1 = ran.choice(species[:initSpeciesListLength-1])
								idsub1 = species.index(sub1)
					
								sub2 = ran.choice(species[:initSpeciesListLength-1])
								idsub2 = species.index(sub2)
								prod = sub1 + sub2
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
						
						if args.creationMethod == 2: rctsToCat = reactionID - 2
						if catalysisID == 0: 
							cats = np.array([[int(catalysisID), int(catalyst), int(rctsToCat), int(0), float(0.5), float(0.25), float(0.5), ran.randint(1,2)]])
							catalysisID += 1
							if args.creationMethod == 2: # IF wim method
								cats = np.vstack([cats,(int(catalysisID), int(catalyst), int(rctsToCat + 1), int(0), float(0.5), float(0.25), float(0.5), ran.randint(1,2))])
								catalysisID += 1
						else: 
							cats = np.vstack([cats,(int(catalysisID), int(catalyst), int(rctsToCat), int(0), float(0.5), float(0.25), float(0.5), ran.randint(1,2))])
							catalysisID += 1
							if args.creationMethod == 2:
								cats = np.vstack([cats,(int(catalysisID), int(catalyst), int(rctsToCat + 1), int(0), float(0.5), float(0.25), float(0.5), ran.randint(1,2))])
								catalysisID += 1
					#print rcts
					#print cats
					#print rcts.shape
					#print cats.shape
					# Create food list
					foodList = range(args.lastFood+1)
					#print "\t\t|- RAF searching step..."
					rafsets = raf.rafComputation(fid_initRafRes, fid_initRafResALL, fid_initRafResLIST, 'tmpDir', prob, averageConn, rcts, cats, foodList, maxlength)
					#print rafsets
					if len(rafsets[2]) > 0: raffound += 1
					#print raffound
					#raw_input("ciao")

					del rcts
					del cats
			
				time2 = time()
				percAcsFounds = raffound/float(iterations)
				timet = time2-time1
				print "\t\t\t %RAF --> ", percAcsFounds, " TIME: ", time2 - time1, " Cleavages: ", nCleavage, " - Condensations: ", nCondensa
				if percAcsFounds >= 0.99: 
					print "\t\t\t Max number of RAFs found"
					increaseYet = False # If RAF is always found so next average conn is assessed
				strToWrite = str(prob) + "\t" + str(maxlength) + "\t" + str(averageConn) + "\t" + str(percAcsFounds) + "\t" + str(timet) + "\n"
				fid_initRafResSUM.write(strToWrite)
				averageConn += 0.1
	# Close text files
	fid_initRafRes.close()
	fid_initRafResLIST.close()
	fid_initRafResALL.close()
	fid_initRafResSUM.close()	
	
		