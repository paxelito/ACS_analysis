#!/usr/bin/env python
# -*- coding: latin-1 -*-
''' MAIN analysis script package for CaRNeSS simulations
	Main file analysis
	
	http://networkx.lanl.gov/reference/readwrite.html
'''

import sys, os # Standard library
import datetime as dt
import linecache as lc
from copy import deepcopy
import glob
from argparse import ArgumentParser
import numpy as np # Scientific library
from numpy import * 

try:
    from pylab import *
except:
    pass
   
from lib.IO import * 
from lib.graph import raf
from lib.graph import network
from lib.dyn import dynamics as dm
import networkx as nx
import matplotlib.pyplot as plt

#æInput parameters definition 
if __name__ == '__main__':
	parser = ArgumentParser(
				description='Graph analysis of the chemistry'
				, epilog='''Graph CHEMISTRY analysis. ''') 
	parser.add_argument('-p', '--strPath', help='Path where files are stored (def: ./)', default='./')	
	parser.add_argument('-f', '--lastFluxID', help='Last ID of the flux species (def: 5)', default='5', type=int)
	parser.add_argument('-m', '--maxDim', help='Max Dimension of the system (def: 6)', default='6', type=int)
	args = parser.parse_args()
	
	# Create absolute paths
	strPath = os.path.abspath(args.strPath)
	tmpDirs = sort(os.listdir(strPath))
	
	# Goes into the simulation folder
	os.chdir(strPath)
	
	# Create stas folders
	ndn = '_0_new_allStatResults'
	newdirAllResults = os.path.join(strPath, ndn)
	if not os.path.isdir(newdirAllResults):
		try: os.mkdir(newdirAllResults)
		except: print "Impossible to create statistic directory", newdirAllResults; sys.exit(1)
	print "\n\n********************\n\n|- Simulation Folder: ", strPath
	# System Type
	_CLOSE_ = 0
	_PROTO_ = 1
	_CSTR_ = 2
	
	fname_initRafRes = os.path.join(newdirAllResults, '0_initRafAnalysis.csv')
	fname_initRafResLIST = os.path.join(newdirAllResults, '0_initRafAnalysisLIST.csv')
	fname_initRafResALL = os.path.join(newdirAllResults, '0_initRafAnalysisALL.csv')
	fid_initRafRes = open(fname_initRafRes, 'w')
	fid_initRafResLIST = open(fname_initRafResLIST, 'w')
	fid_initRafResALL = open(fname_initRafResALL, 'w')
	strToWrite = "Folder\tP\tAC\tM\tRAFsize\tClosureSize\tCatsSize\tuRAF\tSCC\tAutoCat\n"
	fid_initRafRes.write(strToWrite)
	
	for tmpDir in tmpDirs:
		
		totDirName = os.path.join(strPath,tmpDir)
		if os.path.isdir(totDirName):
			# Move to the directory 
			os.chdir(totDirName)
			
			print " \- Results Folder: ", totDirName				
			# Analysis of the initial structures 
			conf = readfiles.readConfFile(totDirName) #ÊConfiguration file upload
			foodList = range(0,args.lastFluxID+1)
							
			# Initial Analysis are turned ON (RAF ANALYSIS)
			print "   |- LOADING init structures..."
			rcts = readfiles.loadAllData(totDirName,'_acsreactions.csv') # reaction file upload
			cats = readfiles.loadAllData(totDirName,'_acscatalysis.csv') #Êcatalysis file upload
			#ÊREAL RAF COMPUTATION (!!!!!! 1 is average connectivity)
			raf, scc, sccg = network.net_analysis_of_static_graphs(fid_initRafRes, fid_initRafResALL, fid_initRafResLIST, tmpDir, conf[9], 1, rcts, cats, foodList, args.maxDim)
			for i in raf: print i
			for i in scc: print i
			for line in nx.generate_adjlist(sccg):
				print(line)
			nx.write_graphml(sccg, os.path.join(newdirAllResults, 'graph.graphml'))
			plt.show()
		
		
					
	fid_initRafRes.close()
	fid_initRafResLIST.close()
	fid_initRafResALL.close()
	
					
				  	
				  	


