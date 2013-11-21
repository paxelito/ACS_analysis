#!/usr/bin/python
# -*- coding: latin-1 -*-
'''Script to perform a RAF (Reflexively Autocatalytic Food Generated) analysis on the initial structures. 
'''
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
      
from lib.IO import readfiles as rf
from lib.dyn import dynamics as dm
from lib.graph import raf 

#ÊInput parameters definition 
if __name__ == '__main__':
	parser = ArgumentParser(
				description='Script to perform a RAF analysis.'
				, epilog='''File with angle trajectories are created. ''') 
	parser.add_argument('-p', '--strPath', help='Path where files are stored', default='./')
	args = parser.parse_args()
	
	# Create absolute paths
	strPath = os.path.abspath(args.strPath)
	# Read the param file and define the environment (Close, Protocell, CSTR)
	_CLOSE_ = 0
	_PROTO_ = 1
	_CSTR_ = 2
	conf = rf.readConfFile(strPath)
	# System type recognition
	if (conf[6] == 0) & (conf[7] > 0): sysType = _PROTO_
	elif (conf[6] > 0) & (conf[7] == 0): sysType = _CSTR_
	elif (conf[6] == 0) & (conf[7] == 0): sysType = _CLOSE_

	# Food list creation (first closure of F)
	closure = dm.generateFluxList(strPath, sysType)
	foodSet = deepcopy(closure)
	# Load reaction and catalysis structures
	rcts = rf.loadAllData(strPath,'_acsreactions.csv')
	cats = rf.loadAllData(strPath,'_acscatalysis.csv')
	# Generate cluser to Food
	print "|- First RAF computation"
	closure = raf.generateClosure(strPath,closure,rcts) # just to params since is made on the init folder. 
	# Check RA condition
	RA = raf.RAcondition(strPath,closure,rcts,cats)
	# Check F condition
	RAF = raf.Fcondition(strPath,closure,RA,rcts)
	RAFlpre = len(RAF)
	# Temporary RAF set
	redRcts = rcts[RAF,:]
	
	print "|- Iterative procedure starts....."
	# If RAF set is not empty the iterative procedure can start
	if len(RAF) > 1:
		RAFlpost = 0
		while (len(RAF) > 0) & (RAFlpre > RAFlpost):
			RAFlpre = len(RAF)
			closure = raf.generateClosure(strPath,foodSet,redRcts)
			#print "CLOSURE ->", closure
			RA = raf.RAcondition(strPath,closure,rcts,cats)
			#print "RA SET ->", RA
			RAF = raf.Fcondition(strPath,closure,RA,rcts)
			#print "RAF set ->", RAF
			redRcts = rcts[RAF,:]
			RAFlpost = len(RAF)
	
	print "|- CLOSURE -> ", closure
	#print "|- RA set -> ",  RA
	print "|- RAF set -> ",  RAF
	
	catalists = raf.findCatforRAF(cats, RAF, closure)
	
	print "|- CATALYSTS -> ",  catalists