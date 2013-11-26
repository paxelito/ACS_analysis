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
   
   
def removeRareRcts(graph, dt, life, nrg, deltat):
	if shape(graph)[0] > 1:
		graph[:,life] = graph[:,life] + deltat
		graph[:,nrg] = graph[:,dt] - graph[:,life]
		graph = graph[graph[:,nrg]>0,:]
		return graph
	else:
		return graph

def fixCondensationReaction(m1, m2, m3, rcts):
	
	if sum((rcts[:,2] == m1) & (rcts[:,3] == m2) & (rcts[:,4] == m3)) == 1:
		#print "- Right RCT"
		#print rcts[((rcts[:,2] == m1) & (rcts[:,3] == m2) & (rcts[:,4] == m3)),:]
		#raw_input("ecco...")
		return m1, m2, m3
	elif sum((rcts[:,2] == m1) & (rcts[:,3] == m3) & (rcts[:,4] == m2)) == 1:
		#print "- REVERSE RCT"
		#print rcts[((rcts[:,2] == m1) & (rcts[:,3] == m3) & (rcts[:,4] == m2)),:]
		#raw_input("ecco...")
		return m1, m3, m2
	else: 
		print r
		print "ERROR!!!!"
		sys.exit(1)
		
	
