#!/usr/bin/env python
# -*- coding: latin-1 -*-

'''Function to analyze angles in time between successive generations. 
'''

import sys, os # Standard library
import datetime as dt
import linecache as lc
import glob
import numpy as np # Scientific library
from numpy import * 
from argparse import ArgumentParser
from lib.visual import graphics as gr
from lib.IO import *

try:
    from pylab import *
except:
    pass

#--------------------------------------------------------------------------------------
# Function to create string zero string vector before graph filename. 
# According to the total number of reactions N zeros will be add before the instant reaction number 
# (e.g. reaction 130 of 10000 the string became '00130')
def zeroBeforeStrNum(tmpl, tmpL):
	strZero = ''
	nZeros = len(str(tmpL)) - len(str(tmpl))
	if nZeros > 0:
		for i in range(0,nZeros): strZero = strZero + '0'
	return strZero

# Input parameters definition 
if __name__ == '__main__':
	parser = ArgumentParser(
				description='ELDA (Emergent Levell Dynamic	 Assessment).'
				, epilog='''Copyright Alessandro Filisetti 2013 ''') 
	parser.add_argument('-p', '--StrPath', help='Simulation outcomes path', default='./')
	parser.add_argument('-m', '--maxLength', help='Max length of the system', default='6', type=int)
	parser.add_argument('-f', '--lastFoodMol', help='Max length of the system', default='5', type=int)
	args = parser.parse_args()
	
	# Convert path in absolute path
	StrPath = os.path.abspath(args.StrPath)
	print "|- Dicrectory -> ", StrPath
	
	resDirPath = os.path.abspath(os.path.join("./", "res"))
	
	newdirAllResults = os.path.abspath(os.path.join(os.curdir, '_0_angles'))
	if not os.path.isdir(newdirAllResults):
		try:
			os.mkdir(newdirAllResults)
		except:
			print "Impossible to create statistic directory", newdirAllResults; sys.exit(1)
	
	print " |- Results Folder-> ", resDirPath
	if os.path.isdir(resDirPath):
		os.chdir(resDirPath)
		param = readfiles.read_sims_conf_file('./acsm2s.conf')
		# Compute delta t of the observed data 
		deltat = param[2] / param[40]
		
		# Find the number of generations
		numberOfGens = len(glob.glob(os.path.join(resDirPath,'timeSpeciesAmount_*')))
		# Creating amount files list
		strAmount = 'timeSpeciesAmount_*'
		amountFiles = sorted(glob.glob(os.path.join(resDirPath,strAmount)))
		print amountFiles
		
		npvector = []
		for gen in range(numberOfGens):
			strZeros = zeroBeforeStrNum(gen+1, numberOfGens)
			fileName = 'timeSpeciesAmount_' + strZeros + str(gen+1) + '_1.csv'
			data = np.loadtxt(open(fileName,"rb"),delimiter="\t")
			npvector.append(data)
			
		# Create a 3D T X N X N data structure where angles will be stored
		cubeangle = np.zeros((npvector[0].shape[0],numberOfGens,numberOfGens))
		# define species limits [nofood, --- , nocomplexes]
		firstspeciesid = args.lastFoodMol + 1
		lastspeciesid = 2**(args.maxLength)-3
		# for each time and for each generation
		for t in range(npvector[0].shape[0]):
			for idgx, gx in enumerate(npvector):
				for idgy, gy in enumerate(npvector):
					vecX = np.array(gx[t,firstspeciesid:lastspeciesid])
					vecY = np.array(gy[t,firstspeciesid:lastspeciesid])
					# Compute coseno
					tmpCos = float(np.dot(vecX,vecY) / (np.linalg.norm(vecX) * np.linalg.norm(vecY)))
					if tmpCos >= 1.0: tmpCos = 1.0
					aseno = np.arccos(tmpCos)  * 360.00 / 6.28
					if aseno < 1e-5: aseno = 0
					cubeangle[t,idgx,idgy] = aseno
		
		# Save angle matrices	
		summarystats = np.zeros((npvector[0].shape[0],2))	
		for t in range(cubeangle.shape[0]):
			print '\t|- Time ', t 
			strZeros = zeroBeforeStrNum(t, cubeangle.shape[0])
			savename = 'angle_' + strZeros + str(t) + '.csv'
			np.savetxt(os.path.join(newdirAllResults,savename), cubeangle[t], fmt='%.4f', delimiter='\t')
			# compute mean and standard deviation
			summarystats[t,0] = np.mean(cubeangle[t])
			summarystats[t,1] = np.std(cubeangle[t])
			savename = 'hist_' + strZeros + str(t) + '.png'
			hi = []
			for x in range(cubeangle[t].shape[0]):
				for y in range(cubeangle[t].shape[1]):
					if cubeangle[t][x,y] > 0: hi.append(cubeangle[t][x,y])
			gr.PlotHist(os.path.join(newdirAllResults,savename), hi, 'Angle', 'Frequency',num_bins=10)
				
		np.savetxt(os.path.join(newdirAllResults,'0_summary.csv'), summarystats, fmt='%.4f', delimiter='\t')
		
		
			
		
					
		
		
			
		
	else:
		print "|- NO RESULTS FOLDER"
		sys.exit(1)
	
	