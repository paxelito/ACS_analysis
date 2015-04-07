#!/usr/bin/python
# -*- coding: latin-1 -*-

import sys, os # Standard library
import glob
from copy import deepcopy
import numpy as np # Scientific library
import networkx as nx
from numpy import * 
from argparse import ArgumentParser
try:
    from pylab import *
except:
    pass
   
def createSimpleGraph(tmpCstr, tmpCats):
	"\t\t|- Cat -> Prod graph creation..."
	g = []
	for id, cat in enumerate(tmpCats):
		if id == 0:
			g = np.array([[int(cat[0]), int(tmpCstr[tmpCstr[:,0] == cat[1],0])]])
		else:
			g = np.vstack([g,(int(cat[0]), int(tmpCstr[tmpCstr[:,0] == cat[1],0]))])
	#print g
	#raw_input("prompt")
	return g

def createNetXGraph(tmpCstr, tmpCats):
	"\t\t|- Cat -> Prod graph creation..."
	Gcatpro = nx.DiGraph()
	#print tmpCstr
	#print tmpCats
	for id, cat in enumerate(tmpCats):
		try:
			if int(tmpCstr[tmpCstr[:,0] == cat[2],1]) == 1:
				Gcatpro.add_weighted_edges_from([(int(cat[1]),int(tmpCstr[tmpCstr[:,0] == cat[2],3]),1)])
				if int(tmpCstr[tmpCstr[:,0] == cat[2],3]) is not int(tmpCstr[tmpCstr[:,0] == cat[2],4]):
					Gcatpro.add_weighted_edges_from([(int(cat[1]),int(tmpCstr[tmpCstr[:,0] == cat[2],4]),1)])
			else:
				Gcatpro.add_weighted_edges_from([(int(cat[1]),int(tmpCstr[tmpCstr[:,0] == cat[2],2]),1)])
			#print Gcatpro.edges()
		except:
			print tmpCats
			print tmpCstr
			print cat
			print tmpCstr[tmpCstr[:,0] == cat[2],1]
			sys.exit(1)
	return Gcatpro

def createNetXGraphForRAF(tmpCstr, tmpClosure, tmpCats):
	"\t\t|- Cat -> Prod graph creation..."
	Gcatpro = nx.DiGraph()
	# Extract catalysts catalysing reactions of the RAF set. 
	for id, cat in enumerate(tmpCats): # For each catalysis
		if cat[1] in tmpClosure: # IF the catalyst is in the closure
			if sum(tmpCstr[:,0] == cat[2]) > 0:# if the are reactions catalyzed by the catalyst
				if int(tmpCstr[tmpCstr[:,0] == cat[2],1]) == 1: 
					Gcatpro.add_weighted_edges_from([(int(cat[1]),int(tmpCstr[tmpCstr[:,0] == cat[2],3]),1)])
					if int(tmpCstr[tmpCstr[:,0] == cat[2],3]) is not int(tmpCstr[tmpCstr[:,0] == cat[2],4]):
						Gcatpro.add_weighted_edges_from([(int(cat[1]),int(tmpCstr[tmpCstr[:,0] == cat[2],4]),1)])
				else:
					Gcatpro.add_weighted_edges_from([(int(cat[1]),int(tmpCstr[tmpCstr[:,0] == cat[2],2]),1)])
	return Gcatpro

def diGraph_netX_stats(tmpDig):
	realSccs = 0
	scc = nx.strongly_connected_components(tmpDig)
	sccsg = nx.strongly_connected_component_subgraphs(tmpDig)
	actualScc = []
	for i in scc: 
		if len(i) > 1: actualScc.append(i)
	sccN= len(actualScc)
	selfLoops = tmpDig.number_of_selfloops()
	selfLoopsEgdes = tmpDig.selfloop_edges()
	realSccs = selfLoops + sccN 
	return actualScc, sccN, selfLoops, selfLoopsEgdes, realSccs, sccsg

def checkMinimalSCCdimension(tmpDig,tmpMinDim):
	findflag = False
	sccSize = []
	scc = nx.strongly_connected_components(tmpDig)
	#print scc
	#raw_input("cioa")
	#sccL = nx.strongly_connected_components_recursive(tmpDig)
	#print scc
	#print sccL
	autocat = tmpDig.nodes_with_selfloops()
	if len(scc) > 0:
		for i in scc: 
			if len(i) >= tmpMinDim:
				if len(i) > 1:
					findflag = True
					sccSize.append(len(i))
				else:
					if i[0] in autocat:
						findflag = True
						sccSize.append(len(i))
				break
			
	return [findflag, sccSize, scc, autocat]

def printSCConFile(tmpSCCL, tmpfolderName, filePrefix):
	filename = '0_' + str(filePrefix) + '_SCCs.txt'
	fname_scc = os.path.join(tmpfolderName, filename)
	fid_fname_scc = open(fname_scc, 'w')
	for sccid, SCC in enumerate(tmpSCCL):
		if len(SCC) >= 1:
			str2w = '\n*** Chemistry ' + str(filePrefix) + ' - SCC number ' + str(sccid+1) + '\n\n'
			fid_fname_scc.write(str2w)
			for scc in SCC: 
				strToWrite = str(scc) + "\t"
				fid_fname_scc.write(strToWrite)
			fid_fname_scc.write("\n")
		




