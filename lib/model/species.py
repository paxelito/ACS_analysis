#!/usr/bin/python
# -*- coding: latin-1 -*-

import sys, os # Standard library
import itertools as it

def createCompleteSpeciesPopulation(M, alphabet):
	species = []
	for i in range(M): species.extend(map(''.join,it.product(alphabet, repeat=i+1)))
	return species

def getTotNumberOfSpeciesFromCompletePop(M):
	N = 2 ** (M + 1) - 2
	return N