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
from ..IO import readfiles as rf
from ..dyn import dynamics as dm
   
def generateClosure(tmpF,rcts):
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

def RAcondition(tmpCL, rcts, cats):
	RAset = []
	for cat in cats:
		if int(cat[1]) in tmpCL: RAset.append(int(cat[2]))
	return RAset

def Fcondition(tmpCL, tmpRA, rcts):
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

# FUNCTION TO FIND RAF IN INIT STRUCTURES 
def rafsearch(conf, rcts, cats, closure):
	
	# Read the param file and define the environment (Close, Protocell, CSTR)
	print "   |- RAF searching step..."
	_CLOSE_ = 0
	_PROTO_ = 1
	_CSTR_ = 2
	# System type recognition
	if (conf[6] == 0) & (conf[7] > 0): sysType = _PROTO_
	elif (conf[6] > 0) & (conf[7] == 0): sysType = _CSTR_
	elif (conf[6] == 0) & (conf[7] == 0): sysType = _CLOSE_

	# Food list creation (first closure of F)
	foodSet = deepcopy(closure)
	# Load reaction and catalysis structures
	# Generate cluser to Food
	# Check RA condition
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
			closure = generateClosure(foodSet,redRcts)
			#print "CLOSURE ->", closure
			RA = RAcondition(closure,rcts,cats)
			#print "RA SET ->", RA
			RAF = Fcondition(closure,RA,rcts)
			#print "RAF set ->", RAF
			redRcts = rcts[RAF,:]
			RAFlpost = len(RAF)
	
	catalists = findCatforRAF(cats, RAF, closure)
	return closure, RA, RAF, catalists
	