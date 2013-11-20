#!/usr/bin/python
# -*- coding: latin-1 -*-

import sys, os # Standard librar
import glob
import numpy as np # Scientific library
from numpy import * 
from argparse import ArgumentParser
try:
    from pylab import *
except:
    pass

# Persona Libraries
from ..IO import readfiles as rf
   
def generateClosure(tmpPath, tmpF,rcts):
	closure = tmpF
	for sngf in closure:
		# Find reduced rct matrix where the species is a reactant 
		reducedR = rcts[((rcts[:,1] == 1)&(rcts[:,2] == sngf))|((rcts[:,1] == 0)&(rcts[:,3] == sngf)),:]
		# Append product from cleavage reaction
		if reducedR.shape[0] > 0:
			for r in reducedR: 
				if int(r[1]) == 1: 
					if int(r[3]) not in closure: closure.append(int(r[3]))
					if int(r[4]) not in closure: closure.append(int(r[4]))
				else:
					if (int(r[3]) in closure) & (int(r[4]) in closure):
						if int(r[2]) not in closure: closure.append(int(r[2]))			
		
	return sorted(closure)

def RAcondition(tmpPath, tmpCL, rcts, cats):
	RAset = []
	for cat in cats:
		if int(cat[1]) in tmpCL: RAset.append(int(cat[2]))
	return RAset

def Fcondition(tmpPath, tmpCL, tmpRA, rcts):
	RAFset = []
	for r in tmpRA:
		if rcts[r,1] == 1:
			if int(rcts[r,2]) in tmpCL: RAFset.append(r)
		else:
			if (int(rcts[r,3]) in tmpCL) & (int(rcts[r,4]) in tmpCL): RAFset.append(r)
	return RAFset

def findCatforRAF(tmpCat, tmpRAF, tmpClosure):
	catalysts = []
	for c in tmpCat:
		if c[2] in tmpRAF:
			if c[1] in tmpClosure:
				catalysts.append(int(c[1]))
	return catalysts