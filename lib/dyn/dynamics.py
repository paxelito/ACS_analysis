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
from ..IO import readfiles as rf
   

   
def generateFluxList(tmpPath, tmpSysType, tmpLastID=None):
	
	tempFood = []
	
	if tmpSysType == 0:
		if tmpLastID: tempFood = range(tmpLastID+1)
		else: print '|- WARNING: No initial flux is generated.'
	elif tmpSysType == 1:
		tempFood = rf.readBufferedID(tmpPath)
	else : 
		tempFood = rf.readCSTRflux(tmpPath)
	
	return tempFood

def rangeFloat(start, step, stop):
	r = start
	while r <= stop:
		yield r
		r += step
		
def fluxAnalysis(tmpDir, resDirPath, strZeros, ngen):
	print "ciao"
		
