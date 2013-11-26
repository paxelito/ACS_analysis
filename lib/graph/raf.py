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
def generateClosure(tmpF,rcts,debug=False):
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
 	if debug:
 		print "1) CLOSURE"
	 	print rcts
	 	print "\t FOOD: ", tmpF
	 	print "\t CLosure: ", sorted(closure)
	return sorted(closure)

# COMPUTE RA CONDITION
def RAcondition(tmpCL, rcts, cats,debug=False):
	RAset = []
	catas = []
	for cat in cats:
		if int(cat[1]) in tmpCL: 
			RAset.append(int(cat[2]))
			catas.append(int(cat[1]))
	if debug:
		print "2) RA"
		print "\t CLosure: ", tmpCL
		print "\t Catalysts selected:", catas
	return RAset

# COMPUTE F CONDITION
def Fcondition(tmpCL, tmpRA, rcts,debug=False):
	RAFset = []
	#print rcts
	#print tmpRA
	try:
		for r in tmpRA:
			if rcts[rcts[:,0]==r,1] == 1:
				if int(rcts[rcts[:,0]==r,2]) in tmpCL: RAFset.append(r)
			else:
				if (int(rcts[rcts[:,0]==r,3]) in tmpCL) & (int(rcts[rcts[:,0]==r,4]) in tmpCL): RAFset.append(r)
		if debug:
			print "3) RAF"
			print "\t RAF: ", RAFset
			print "\t Closure: ", tmpCL
		return RAFset
	except: 
		print rcts
		print tmpRA
		sys.exit(1)
		

# FIND CATALYSTS OF THE RAF SET
def findCatforRAF(tmpCat, tmpRAF, tmpClosure,debug=False):
	catalysts = []
	for c in tmpCat:
		if c[2] in tmpRAF:
			if c[1] in tmpClosure:
				catalysts.append(int(c[1]))
	return catalysts

# FUNCTION TO FIND RAF IN INIT STRUCTURES 
def rafsearch(rcts, cats, closure,debug=False):
	
	if rcts.shape[0] > 0:
		
		# Food list creation (first closure of F)
		foodSet = deepcopy(closure)
		closure = generateClosure(closure,rcts,debug)
		RA = RAcondition(closure,rcts,cats,debug)
		# Check F condition
		RAF = Fcondition(closure,RA,rcts,debug)
		RAFlpre = len(RAF)
		# Temporary RAF set
		redRcts = rcts[np.any(rcts[:, 0] == np.expand_dims(RAF,1), 0), :]
		if debug:
			print RAF
			print redRcts
		
		# If RAF set is not empty the iterative procedure can start
		if len(RAF) > 1:
			RAFlpost = 0
			while (len(RAF) > 0) & (RAFlpre > RAFlpost):
				RAFlpre = len(RAF)
				foodCopy = deepcopy(foodSet)
				closure = generateClosure(foodCopy,redRcts,debug)
				RA = RAcondition(closure,rcts,cats,debug)
				RAF = Fcondition(closure,RA,rcts,debug)
				redRcts = rcts[np.any(rcts[:, 0] == np.expand_dims(RAF,1), 0), :]
				RAFlpost = len(RAF)
		
		catalists = findCatforRAF(cats, RAF, closure)
		if debug:
			print "4) FINAL"
			print "\t Closure -> ", closure
			print "\t cats -> ", catalists
			print "\t RAF -> ", RAF
			print "\t matx dimension -> ", redRcts.shape
			print "\n************************************\n"
			raw_input("closure stop!!!") 
		return closure, RA, RAF, catalists, len(list(set(RAF)))
	else:
		return [], [], [], [], 0

# BRIDGE FUNCTION TO DETECT RAFs in INITIAL SETS
def rafComputation(fid_initRafRes, fid_initRafResALL, fid_initRafResLIST, tmpDir, rctProb, avgCon, rcts, cats, foodList, maxDim,debug=False):
	rafset = rafsearch(rcts, cats, foodList,debug) # RAF search 
	ErctP = "%.4g" % rctProb
	strToWrite = tmpDir + "\t" + str(ErctP) + "\t" + str(avgCon) + "\t" + str(maxDim) + "\t" + str(len(rafset[2])) + "\t" + str(len(rafset[0])) + "\t" + str(len(rafset[3])) + "\t" + str(rafset[4]) + "\n"
	fid_initRafRes.write(strToWrite)
	writefiles.write_init_raf_list(fid_initRafResLIST, rafset, tmpDir)
	writefiles.write_init_raf_all(fid_initRafResALL, rafset, tmpDir, rcts, cats)
	return rafset

# BRIDGE FUNCTION TO DETECT RAFs in DYNAMICS
def rafDynamicComputation(fid_dynRafRes, tmpTime, rcts, cats, foodList, growth=False, rctsALL=None, catsALL=None, completeRCTS=None,debug=False):
	#print rcts
	#print cats
	rafset = rafsearch(rcts, cats, foodList,debug) # RAF search
	if growth == True: rafsetALL = rafsearch(rctsALL, catsALL, foodList) # RAF search
	strRAF = '' 
	# If RAF analysis is made in dynamical temporary structures a trnaslation in real net must be done
	if completeRCTS != None: convRAF = findRAFrcts(rafset[2],rcts,completeRCTS)
	else: convRAF = rafset[2]
	if len(convRAF) > 0: 		
		for x in convRAF: strRAF += str(x) + '\t'	
	if growth == False: strToWrite = str(tmpTime) + "\t" + str(len(rafset[0])) + "\t" + str(rafset[4]) + "\t" + strRAF + "\n"
	else: strToWrite = str(tmpTime) + "\t" + str(len(rafset[0])) + "\t" + str(rafset[4]) + str(len(rafsetALL[0])) + "\t" + str(rafsetALL[4]) + "\t" + strRAF + "\n"
	fid_dynRafRes.write(strToWrite)
	return rafset

def findRAFrcts(RAF, rcts, actrcts):
	rafset = []
	for i in RAF:
		position = ((actrcts[:,1] == rcts[rcts[:,0]==i,1]) & (actrcts[:,2] == rcts[rcts[:,0]==i,2]) & (actrcts[:,3] == rcts[rcts[:,0]==i,3]))
		rafset.append(int(actrcts[position,0]))
	return list(set(rafset))

	
	