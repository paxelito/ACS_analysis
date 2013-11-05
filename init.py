#!/usr/bin/env python
# -*- coding: latin-1 -*-
'''Script to order the analysis of the divergences in time. 
'''
import sys, os # Standard librar
import glob
import numpy as np # Scientific library
import itertools as it
import random
from numpy import * 
from argparse import ArgumentParser
try:
    from pylab import *
except:
    pass

#ÊInput parameters definition 
if __name__ == '__main__':
	parser = ArgumentParser(
				description='This script re-arrange results in a more friendly way from the angle analysis in time.'
				, epilog='''File with angle trajectories are created. ''') 
	parser.add_argument('-s', '--strPath', help='Path where files are stored', default='./')
	parser.add_argument('-m', '--maxDim', help='Max Dimension of the system', default='6', type=int)
	parser.add_argument('-p', '--rctProb', help='Reaction Probability', default='0', type=float)
	parser.add_argument('-c', '--noCat', help='Non catalysis max length', default='2', type=float)
	parser.add_argument('-r', '--randomSeed', help='random seed', type=int, default=None)
	args = parser.parse_args()
	
	print "\nSystem Initialization..."
	np.random.seed(args.randomSeed)
	# Create all the the species starting from the alphabet
	alphabet = ['A', 'B']
	species = []
	for i in range(args.maxDim): species.extend(map(''.join,it.product(alphabet, repeat=i+1)))
	
	# compute number of cleavage
	totSpecies = 2 ** (args.maxDim + 1) - 2
	print '|- Total Number of Species: ', totSpecies
	# Compute overall conceivable number of reactions
	totCleavage = sum(map(lambda x: len(x)-1,species))
	totCond = totSpecies ** 2
	totRcts = totCleavage + totCond
	
	# If the reaction probability is 0, it is set to the critical value of 1 reaction per species on average
	if args.rctProb == 0: prob = 1 / float(totSpecies) 
	else: prob = args.rctProb
	
	print '|- Total Number of Reactions: ', totRcts, ' <- Condensations: ', totCond, ' - Cleavage: ', totCleavage
	# Compute reactions to catalyse according to the probability
	rctToCat = int(round(totRcts * totSpecies * prob))
	print '|- According to probability ', prob, ' and the number of species (', totSpecies,'), ', rctToCat, ' catalysis will be created. '
	
	for i in range(rctToCat):
		# Select if condensation of cleavage according to the total number of reactions
		rnd = np.random.randint(totRcts)
		rctType = 0
		if rnd < totCleavage: rctType = 1
		
		# If cleavage
		if rctType == 1:
			# Select species to cleave
			print "To be continued..."
		
		print rnd, " ", totRcts, " ", totCleavage, " ", rctType
	
	print '\n'
	
	
		