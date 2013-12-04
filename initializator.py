	#!/usr/bin/env python
# -*- coding: latin-1 -*-
'''Script to order the analysis of the divergences in time. 
'''
import sys, os # Standard librar
import glob
import numpy as np # Scientific library
import itertools as it
import random as ran
from time import time
from numpy import * 
from argparse import ArgumentParser
try:
    from pylab import *
except:
    pass
   
from lib.graph import raf
from lib.graph import scc
from lib.graph import network
from lib.dyn import dynamics as dn
from lib.model import species as sp
from lib.model import reactions
from lib.IO import *


#æInput parameters definition 
if __name__ == '__main__':
	parser = ArgumentParser(
				description='This script re-arrange results in a more friendly way from the angle analysis in time.'
				, epilog='''File with angle trajectories are created. ''') 
	parser.add_argument('-a', '--sysType', help='System Architecture [1:CLOSE, 2:PROTO, 3:CSTR], deafult: 1', default='1')
	parser.add_argument('-k', '--creationMethod', help='Network creation method (1: Filisetti, 2: Wim, 3: WimNoRevs, DEF: 1)', default='1', type=int)
	parser.add_argument('-f', '--lastFood', type=int, help='Last food species ID (deafult: 5)', default='2')
	parser.add_argument('-n', '--noCat', help='Non catalytic max length (default: 2)', default='2', type=int)
	parser.add_argument('-s', '--initSet', type=int, help='Dimension of the initial set (Default: 4)', default='4')
	parser.add_argument('-m', '--maxDim', help='Max Dimension of the systems (Default: 6)', default='6', type=int)
	parser.add_argument('-o', '--strOut', help='Path for output file storing (Default: ./)', default='./')
	parser.add_argument('-I', '--conf', help='Configuration File (Default: ./acsm2s.conf)', default='./acsm2s.conf')
	parser.add_argument('-i', '--iteration', help='How many network instances to created (Default: 4)', default='4', type=int)
 	parser.add_argument('-v', '--avgCon', help='Catalysis level (deafult: 1), i.e. average reactions catalyzed per species', type=int, default='1')
	parser.add_argument('-c', '--rctRatio', help='Ratio between cleavages and condensations (default: 0, it means that the actual ratio is used)', default='0', type=float)
	parser.add_argument('-C', '--core', help='Number of core on which simulations are distributed', default='2', type=int)	
	parser.add_argument('-F', '--folderName', help='Simulation Folder Name (Deafault: SIMS)', default='SIMS')	
	parser.add_argument('-r', '--randomSeed', help='random seed', type=int, default=None)
	args = parser.parse_args()
	
	# SIMULATION FOLDER CREATION
	folderName = os.path.join(args.strOut, args.folderName)
	if not os.path.isdir(folderName):
		try:
			os.mkdir(folderName)
		except:
			print "Impossible to create statistic directory", folderName; sys.exit(1)
			
	# FOR EACH NETWORK
	fid_run = []
	for sngCore in range(args.core):
		zeroBeforeName =  readfiles.zeroBeforeStrNum(sngCore+1, args.core)
		runFileName = zeroBeforeName + str(sngCore+1) + '_simulation.sh'
		fname_run = os.path.join(folderName, runFileName)
		fid_run.append(open(fname_run, 'w'))
		
	# Read Conf File
	parameters = readfiles.read_sims_conf_file(args.conf)
	
	fidid = 0 # Core on which the chemistry will run
	for singlechem in range(args.iteration):
		# CREATE NETWORK FOLDER
		zeroBeforeName =  readfiles.zeroBeforeStrNum(singlechem+1, args.iteration)
		chemFolder = 'CH' + str(zeroBeforeName) + str(singlechem+1)
		chemFolderPath = os.path.join(folderName, chemFolder)
		if not os.path.isdir(chemFolderPath):
			try:
				os.mkdir(chemFolderPath)
			except:
				print "Impossible to create statistic directory", chemFolderPath; sys.exit(1)	
		
		resFolder = os.path.join(chemFolderPath, 'res')
		if not os.path.isdir(resFolder):
			try:
				os.mkdir(resFolder)
			except:
				print "Impossible to create statistic directory", resFolder; sys.exit(1)
		
		# -----------------------------
		# ARTIFICIAL CHEMISTRY CREATION
		# -----------------------------
		# Save config file
		writefiles.write_acsms_file(chemFolderPath,*parameters)
		writefiles.write_and_create_std_nrgFile(chemFolderPath)
		
				
		
		# UPDATE SIMULATION LUNCHER
		str2w = "echo 'Simulation " + chemFolderPath + "'\nnice ./carness ./" + chemFolderPath + "/ " + chemFolderPath + "/res/ " + chemFolderPath + "/ > 0_log_" + chemFolderPath + ".log\n"
		
		fid_run[fidid].write(str2w)
		fidid += 1
		if fidid == args.core: fidid = 0 
	
	# Close fid runs		
	map(close(),fid_run)