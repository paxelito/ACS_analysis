	#!/usr/bin/env python
# -*- coding: latin-1 -*-
'''Script to order the analysis of the divergences in time. 
'''
import sys, os # Standard librar
import glob
import numpy as np # Scientific library
import itertools as it
import random as ran
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
	parser.add_argument('-a', '--sysType', help='System Architecture [1:CLOSE, 2:PROTO, 3:CSTR]', default='1')
	parser.add_argument('-f', '--lastFood', type=int, help='Last food species ID', default='29')
	parser.add_argument('-s', '--strPath', help='Path where files are stored', default='./')
	parser.add_argument('-o', '--strOut', help='Path for output file storing', default='./')
	parser.add_argument('-m', '--maxDim', help='Max Dimension of the systems (Default: 6)', default='6', type=int)
	parser.add_argument('-i', '--iteration', help='How many network instances per dimension are created (Default: 10)', default='10', type=int)
	parser.add_argument('-p', '--rctProb', help='Reaction Probability, if 0 args.avgCon will be used', default='0', type=float)
	#parser.add_argument('-v', '--avgCon', help='Tuple containing the average connectivity range', default=[0.5, 0.05, 1], type=tuple) nargs=2, metavar=('bar', 'baz')
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
	strToWrite = "P\tM\tAC\tRAF%\n"
	fid_initRafResSUM.write(strToWrite)
	
	# Create food list
	foodList = range(args.lastFood+1)
	
	for maxlength in range(3,args.maxDim+1): # For each dimension
		avgCon = dn.rangeFloat(float(args.avgCon[0]), float(args.avgCon[1]), float(args.avgCon[2]))
		if args.rctProb == 0: # if reaction probability is 0
			for averageConn in avgCon: # For each average connectivity
				raffound = 0
				for instanceID, instance in enumerate(range(args.iteration)): # For each instance of network
					# Create all the the species starting from the alphabet
					alphabet = ['A', 'B']
					species = []
					for i in range(maxlength): species.extend(map(''.join,it.product(alphabet, repeat=i+1)))
					
					# compute number of cleavage
					totSpecies = 2 ** (maxlength + 1) - 2
					# Compute overall conceivable number of reactions
					totCleavage = sum(map(lambda x: len(x)-1,species))
					totCond = totSpecies ** 2
					totRcts = totCleavage + totCond
					
					# If the reaction probability is 0, it is set to the critical value of 1 reaction per species on average
					if args.rctProb == 0: prob = (1 / float(totRcts)) * averageConn
					else: prob = args.rctProb
					
					if instanceID == 0: print '|- Max Length: ', maxlength, ' - AVG CON: ', averageConn, ' - Species: ', totSpecies, '- Reactions: ', totRcts, ' <- Condensations: ', totCond, ' - Cleavage: ', totCleavage
					# Compute reactions to catalyse according to the probability
					rctToCat = int(round(totRcts * totSpecies * prob))
					if instanceID == 0: print '\t\t|- According to probability ', prob, ' and the number of species (', totSpecies,'), ', rctToCat, ' catalysis will be created. '
					
					nCleavage = 0
					nCondensa = 0
					initSpeciesListLength = len(species)
					
					#return (ngen,nsim,totTimes,nrgType,totalRcts,nrgConc,influx_rate,maxLOut,lfdsID,rctProb)
					conf = (1,1,2000,0,200000,0,0,2,args.lastFood,prob)	
				
					for i in range(rctToCat):
						# Select if condensation of cleavage according to the total number of reactions
						#if rctToCat > 1000: 
						#	if i % 1000 == 0: print "\t\t\t|- Reaction ", i, " - species list L: ", len(species)
						rctType = 1
						if args.rctRatio > 0: 
							if ran.random() > args.rctRatio: rctType = 0
						else:
							if (totCleavage / (float(totCleavage) + totCond)) <= ran.random(): rctType = 0
						
						# CREATE REACTION 
						# If cleavage
						if (rctType == 1) & (nCleavage <= totCleavage):
							# Select species to cleave
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
							
							# Reaction Structure Creation
							if i == 0: rcts = np.array([[int(i), int(rctType), species.index(molToCleav), tmp1id, tmp2id, int(0), int(0), int(0)]])
							else: rcts = np.vstack([rcts,(int(i), int(rctType), species.index(molToCleav), tmp1id, tmp2id, int(0), int(0), int(0))])		
							nCleavage += 1
						else: # condensation
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
							# Reaction Structure Creation
							if i == 0: rcts = np.array([[int(i), int(rctType), tmpprodid, idsub1, idsub2, int(0), int(0), int(0)]])
							else: rcts = np.vstack([rcts,(int(i), int(rctType), tmpprodid, idsub1, idsub2, int(0), int(0), int(0))])	
							nCondensa += 1
						# A CATALYST IS RANDOMLY ASSIGNED FROM THE SPECIES LIST
						catalyst = -1
						catFound = False
						while not catFound: 
							catalyst = species.index(ran.choice(species[len(alphabet):]))
							if len(species[catalyst]) > args.noCat: catFound = True
						if i == 0: cats = np.array([[int(i), int(catalyst), int(i), int(0), float(0.5), float(0.25), float(0.5), ran.randint(1,2)]])
						else: cats = np.vstack([cats,(int(i), int(catalyst), int(i), int(0), float(0.5), float(0.25), float(0.5), ran.randint(1,2))])
			
					#print "\t\t|- RAF searching step..."
					rafsets = raf.rafComputation(fid_initRafRes, fid_initRafResALL, fid_initRafResLIST, 'tmpDir', prob, averageConn, rcts, cats, foodList, maxlength)
					if len(rafsets[2]) > 0: raffound += 1
					del rcts
					del cats
	
				strToWrite = str(prob) + "\t" + str(maxlength) + "\t" + str(averageConn) + "\t" + str(raffound/float(args.iteration)) + "\n"
				fid_initRafResSUM.write(strToWrite)
	# Close text files
	fid_initRafRes.close()
	fid_initRafResLIST.close()
	fid_initRafResALL.close()
	fid_initRafResSUM.close()	
	
		