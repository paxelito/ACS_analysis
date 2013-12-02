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
   
from ..graph import raf
from ..graph import scc
from ..dyn import dynamics as dn
from ..IO import writefiles
   
   
def removeRareRcts(graph, dt, life, nrg, deltat):
	if shape(graph)[0] > 1:
		graph[:,life] = graph[:,life] + deltat
		graph[:,nrg] = graph[:,dt] - graph[:,life]
		graph = graph[graph[:,nrg]>0,:]
		return graph
	else:
		return graph

def fixCondensationReaction(m1, m2, m3, rcts):
	
	if sum((rcts[:,2] == m1) & (rcts[:,3] == m2) & (rcts[:,4] == m3)) > 1:
		#print "- Right RCT"
		#print rcts[((rcts[:,2] == m1) & (rcts[:,3] == m2) & (rcts[:,4] == m3)),:]
		#raw_input("ecco...")
		return m1, m2, m3
	elif sum((rcts[:,2] == m1) & (rcts[:,3] == m3) & (rcts[:,4] == m2)) > 1:
		#print "- REVERSE RCT"
		#print rcts[((rcts[:,2] == m1) & (rcts[:,3] == m3) & (rcts[:,4] == m2)),:]
		#raw_input("ecco...")
		return m1, m3, m2
	else: 
		print m1, m2, m3
		print "ERROR!!!!"
		sys.exit(1)

# BRIDGE FUNCTION TO DETECT RAFs in INITIAL SETS
def net_analysis_of_static_graphs(fid_initRafRes, fid_initRafResALL, fid_initRafResLIST, tmpDir, rctProb, avgCon, rcts, cats, foodList, maxDim,debug=False):
	rafset = raf.rafsearch(rcts, cats, foodList,debug) # RAF search 
	stdgraph = scc.createNetXGraph(rcts,cats)
	sccsets = scc.diGraph_netX_stats(stdgraph)
	ErctP = "%.4g" % rctProb
	strToWrite = tmpDir + "\t" + str(ErctP) + "\t" + str(avgCon) + "\t" + str(maxDim) + "\t" + \
				 str(len(rafset[2])) + "\t" + str(len(rafset[0])) + "\t" + str(len(rafset[3])) + \
				 "\t" + str(rafset[4]) + "\t" + str(sccsets[4]) + "\t" + str(sccsets[2]) + "\n"
	fid_initRafRes.write(strToWrite)
	writefiles.write_init_raf_list(fid_initRafResLIST, rafset, tmpDir)
	writefiles.write_init_raf_all(fid_initRafResALL, rafset, tmpDir, rcts, cats)
	return rafset, sccsets

# BRIDGE FUNCTION TO DETECT RAFs in DYNAMICS
def net_analysis_of_dynamic_graphs(fid_dynRafRes, tmpTime, rcts, cats, foodList, growth=False, rctsALL=None, catsALL=None, completeRCTS=None,debug=False):
	#print rcts
	#print cats
	rafset = raf.rafsearch(rcts, cats, foodList,debug) # RAF search
	stdgraph = scc.createNetXGraph(rcts,cats) # netX graph creation
	sccsets = scc.diGraph_netX_stats(stdgraph) # SCC analysis
	if growth == True: 
		rafsetALL = raf.rafsearch(rctsALL, catsALL, foodList) # RAF search
		stdgraphALL = scc.createNetXGraph(rctsALL,catsALL) # netX graph creation
		sccsetsALL = scc.diGraph_netX_stats(stdgraphALL) # SCC analysis
		
	strRAF = '' 
	# If RAF analysis is made in dynamical temporary structures a trnaslation in real net must be done
	 
	if completeRCTS != None: convRAF = raf.findRAFrcts(rafset[2],rcts,completeRCTS)
	else: convRAF = rafset[2]
	if len(convRAF) > 0: 		
		for x in convRAF: strRAF += str(x) + '\t'	
	if growth == False: strToWrite = str(tmpTime) + "\t" + str(len(rafset[0])) + "\t" + str(rafset[4]) + "\t" + str(sccsets[4]) + "\t" + str(sccsets[2]) + "\t" + strRAF + "\n"
	else: strToWrite = str(tmpTime) + "\t" + str(len(rafset[0])) + "\t" + str(rafset[4]) + str(len(rafsetALL[0])) + "\t" + str(rafsetALL[4]) + "\t" + str(sccsets[4]) + "\t" + str(sccsets[2]) + "\t" + strRAF + "\n"
	fid_dynRafRes.write(strToWrite)
	return rafset
		
	
