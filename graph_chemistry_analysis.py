#!/usr/bin/env python
# -*- coding: latin-1 -*-
''' 
	This python tool evaluates a particular chemistry findind RAF, SCC and 
	saving the multigraph bipartite network and the catalyst-product network
	
	NETWORKX formats :: http://networkx.lanl.gov/reference/readwrite.html
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
from lib.visual import graphics as grf
import networkx as nx
import matplotlib.pyplot as plt
from networkx.algorithms import bipartite
from lib.graph import scc

#�Input parameters definition 
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
			
			print " \- Results Folder: {0}".format(totDirName)				
			# Analysis of the initial structures 
			conf = readfiles.readConfFile(totDirName) #�Configuration file upload
			foodList = range(0,args.lastFluxID+1)
							
			# Initial Analysis are turned ON (RAF ANALYSIS)
			print "   |- LOADING init structures..."
			rcts = readfiles.loadAllData(totDirName,'_acsreactions.csv') # reaction file upload
			cats = readfiles.loadAllData(totDirName,'_acscatalysis.csv') #�catalysis file upload

			raf, _, sccg = network.net_analysis_of_static_graphs(fid_initRafRes, fid_initRafResALL, fid_initRafResLIST, tmpDir, conf[9], 1, rcts, cats, foodList, args.maxDim)

			grf.plotBipartiteGraph(rcts, cats, newdirAllResults, "0_graph_bipartite.graphml", "completebipartitegraph.png")

			if len(raf[2]) > 0:
				# Filter graf network
				rafcats = cats[np.in1d(cats[:,1], raf[3])]
				rafrcts = rcts[np.in1d(rcts[:,0], raf[2])]
				grf.plotBipartiteGraph(rafrcts, rafcats, newdirAllResults, "0_graph_RAF_bipartite.net", "bipartiteRAFgraph.png", True)
				grf.plotGraph(sccg, newdirAllResults, "0_graph_catalyst_product.net", "chemistrygraph.png", True)

					
	fid_initRafRes.close()
	fid_initRafResLIST.close()
	fid_initRafResALL.close()
	
	