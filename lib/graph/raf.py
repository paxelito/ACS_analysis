#!/usr/bin/python
# -*- coding: latin-1 -*-

import sys, os # Standard librar
import glob
from copy import deepcopy
import numpy as np # Scientific library
from numpy import * 
from argparse import ArgumentParser
try:
    from pylab import *
except:
    pass

# Persona Libraries
from ..IO import *
from ..dyn import dynamics as dm
   
# GENERATE CLOSURE ACCORDING TO THE FOODSET AND THE SPECIFIC SUBSET OF REACTIONS
def generateClosure(tmpF,rcts):
	closure = tmpF
	for sngf in closure:
		#print "|- MOL IN CLO: ", sngf
		# Find reduced rct matrix where the species is a reactant 
		reducedR = rcts[((rcts[:,1] == 1)&(rcts[:,2] == sngf))|(((rcts[:,1] == 0)&(rcts[:,3] == sngf))|((rcts[:,1] == 0)&(rcts[:,4] == sngf))),:]
		# Append product from cleavage reaction
		if reducedR.shape[0] > 0:
			for r in reducedR: 
				#print "\t|- rct: ",  r
				if int(r[1]) == 1: 
					if int(r[3]) not in closure: closure.append(int(r[3]))
					if int(r[4]) not in closure: closure.append(int(r[4]))
				else:
					if (int(r[3]) in closure) & (int(r[4]) in closure):
						if int(r[2]) not in closure: closure.append(int(r[2]))	
				#print "\t|- closure: ",   closure
			#raw_input("...")	

	return sorted(closure)

# COMPUTE RA CONDITION
def RAcondition(tmpCL, rcts, cats):
	RAset = []
	for cat in cats:
		if int(cat[1]) in tmpCL: RAset.append(int(cat[2]))
	return RAset

# COMPUTE F CONDITION
def Fcondition(tmpCL, tmpRA, rcts):
	RAFset = []
	for r in tmpRA:
		if rcts[r,1] == 1:
			if int(rcts[r,2]) in tmpCL: RAFset.append(r)
		else:
			if (int(rcts[r,3]) in tmpCL) & (int(rcts[r,4]) in tmpCL): RAFset.append(r)
	return RAFset

# FIND CATALYSTS OF THE RAF SET
def findCatforRAF(tmpCat, tmpRAF, tmpClosure):
	catalysts = []
	for c in tmpCat:
		if c[2] in tmpRAF:
			if c[1] in tmpClosure:
				catalysts.append(int(c[1]))
	return catalysts

# FUNCTION TO FIND RAF IN INIT STRUCTURES 
def rafsearch(rcts, cats, closure):
	
	# Food list creation (first closure of F)
	foodSet = deepcopy(closure)
	closure = generateClosure(closure,rcts)
	RA = RAcondition(closure,rcts,cats)
	# Check F condition
	RAF = Fcondition(closure,RA,rcts)
	RAFlpre = len(RAF)
	# Temporary RAF set
	redRcts = rcts[RAF,:]
	
	# If RAF set is not empty the iterative procedure can start
	if len(RAF) > 1:
		RAFlpost = 0
		while (len(RAF) > 0) & (RAFlpre > RAFlpost):
			RAFlpre = len(RAF)
			foodCopy = deepcopy(foodSet)
			closure = generateClosure(foodCopy,redRcts)
			RA = RAcondition(closure,rcts,cats)
			RAF = Fcondition(closure,RA,rcts)
			redRcts = rcts[RAF,:]
			RAFlpost = len(RAF)
	
	catalists = findCatforRAF(cats, RAF, closure)
	return closure, RA, RAF, catalists, len(list(set(RAF)))

# BRIDGE FUNCTION TO DETECT RAFs
def rafComputation(fid_initRafRes, fid_initRafResALL, fid_initRafResLIST, tmpDir, rctProb, avgCon, rcts, cats, foodList, maxDim):
	rafset = rafsearch(rcts, cats, foodList) # RAF search 
	ErctP = "%.4g" % rctProb
	strToWrite = tmpDir + "\t" + str(ErctP) + "\t" + str(avgCon) + "\t" + str(maxDim) + "\t" + str(len(rafset[2])) + "\t" + str(len(rafset[0])) + "\t" + str(len(rafset[3])) + "\t" + str(rafset[4]) + "\n"
	fid_initRafRes.write(strToWrite)
	writefiles.write_init_raf_list(fid_initRafResLIST, rafset, tmpDir)
	writefiles.write_init_raf_all(fid_initRafResALL, rafset, tmpDir, rcts, cats)
	return rafset

	
	