#!/usr/bin/python
# -*- coding: latin-1 -*-
'''This python program assess different network structures in term of topological (NON DYNAMICAL) RAF and SCC presence according 
	to different structural parameters. 
'''
import sys, os # Standard library
import glob
from argparse import ArgumentParser
import numpy as np # Scientific library
from numpy import * 
from lib.visual import graphics as gr

try:
    from pylab import *
except:
    pass


if __name__ == '__main__':
	parser = ArgumentParser(
				description='This script create graphs with regard to the network topological assessment.'
				, epilog='''Random Catalytic Reaction Networks Topological Assessment graphic view''') 
	parser.add_argument('-p', '--strPath', help='Path where files are stored (def: ./)', default='./')
	parser.add_argument('-d', '--minDim', help='Minimal Dimension Analysis (def: 5)', default='5', type=int)
	parser.add_argument('-D', '--maxDim', help='Maximal Dimension Analysis (def: 8)', default='8', type=int)
	parser.add_argument('-f', '--imgFormat', help='Image Format', default='eps')
	
	args = parser.parse_args()
	
	# Create absolute paths
	strPath = os.path.abspath(args.strPath)
	tmpDirs = sort(os.listdir(strPath))
	
	for tmpdir in tmpDirs: # For each folder
		if os.path.isdir(os.path.join(strPath, tmpdir)): # If it is an actual folder
			print "|- Graph in ", tmpdir
			print "\t|- Loading row data from summary file..."
			# move into the specific directory
			os.chdir(os.path.join(strPath,tmpdir))
			# load data
			dat = np.loadtxt(open('0_initRafAnalysisSUM.csv'), delimiter='\t', skiprows=1)
			# Create plot
			print "\t|- Create RAF graph..."
			tmpfilename = '1_raf_' + tmpdir
			gr.plot_np_txt_selcols_and_rows_allvarinrows(tmpfilename,dat,(2,3),range(args.minDim,args.maxDim+1),1,'Level of catalysis','RAFs',\
														 range(args.minDim,args.maxDim+1),(0,1),None,tmpimgformat=args.imgFormat)
			print "\t|- Create SCC graph..."
			tmpfilename = '2_scc_' + tmpdir
			gr.plot_np_txt_selcols_and_rows_allvarinrows(tmpfilename,dat,(2,4),range(args.minDim,args.maxDim+1),1,'Level of catalysis','SCCs',\
														 range(args.minDim,args.maxDim+1),(0,1),None,tmpimgformat=args.imgFormat)
			
			print "\t|- Create SCC in RAF graph..."
			tmpfilename = '3_sccraf_' + tmpdir
			gr.plot_np_txt_selcols_and_rows_allvarinrows(tmpfilename,dat,(2,5),range(args.minDim,args.maxDim+1),1,'Level of catalysis','SCCs in RAFs',\
														 range(args.minDim,args.maxDim+1),(0,1),None,tmpimgformat=args.imgFormat)





			
		