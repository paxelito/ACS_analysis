#!/usr/bin/python
# -*- coding: latin-1 -*-

import sys, os # Standard library
import random as ran

def getNumOfCleavages(tmpSpecies):
	return sum(map(lambda x: len(x)-1,tmpSpecies))

def getNumOfCondensations(N):
	return N ** 2

def createRandomCleavage(tmpSpecies, alphabet, tmpInitLMax):
	M1 = ran.choice(tmpSpecies[len(alphabet):tmpInitLMax-1])
	cutPt = ran.randint(1,len(M1)-1) 
	M2 = M1[0:cutPt]
	
	try: 
		M2_ID = tmpSpecies.index(M2)
		M2new = True
	except: 
		M2_ID = len(tmpSpecies)+1
		M2new = False
		
	M3 = M1[cutPt:len(M1)]
	try: 
		M3_ID = tmpSpecies.index(M3)
	except: 
		if M2new == True: M3_ID = len(tmpSpecies) + 1   
		else: M3_ID = len(tmpSpecies) + 2
		
	return M1, M2, M3, M2_ID, M3_ID, M2new

def createRandomCondensation(tmpSpecies, tmpInitLMax):
	sub1 = ran.choice(tmpSpecies[:tmpInitLMax-1])
	idsub1 = tmpSpecies.index(sub1)

	sub2 = ran.choice(tmpSpecies[:tmpInitLMax-1])
	idsub2 = tmpSpecies.index(sub2)
	prod = sub1 + sub2
	
	return sub1, sub2, idsub1, idsub2, prod 