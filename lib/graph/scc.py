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
		Gcatpro.add_weighted_edges_from([(int(cat[1]),int(tmpCstr[tmpCstr[:,0] == cat[2],0]),1)])
	return Gcatpro

def diGraph_netX_stats(tmpDig):
	realSccs = 0
	scc = nx.strongly_connected_components(tmpDig)
	actualScc = []
	for i in scc: 
		if len(i) > 1: actualScc.append(i)
	sccN= len(actualScc)
	selfLoops = tmpDig.number_of_selfloops()
	selfLoopsEgdes = tmpDig.selfloop_edges()
	realSccs = selfLoops + sccN 
	return actualScc, sccN, selfLoops, selfLoopsEgdes, realSccs



