#!/usr/bin/python
# -*- coding: latin-1 -*-

import sys, os # Standard librar
import glob
import numpy as np # Scientific library
import random as ran
import cPickle as pickle
from numpy import * 
from argparse import ArgumentParser
try:
    from pylab import *
except:
    pass
   
def saveRandomSeed(tmpPath):
	'''Function to save the random seed'''
	
	with open(os.path.join(tmpPath,"rndstate.dat"), 'wb') as f:
		pickle.dump(ran.getstate(), f)

def write_init_raf_list(fid, rafinfo, folder):
	fid.write(folder + "\nCLOSURE -> ")
	map(lambda x: fid.write((str(int(x)) + ' ')), rafinfo[0])
	fid.write("\nRAF -> ")
	map(lambda x: fid.write((str(int(x)) + ' ')), rafinfo[2])
	fid.write("\nCATALYSTS -> ")
	map(lambda x: fid.write((str(int(x)) + ' ')), rafinfo[3])
	fid.write("\n\n")
	
def write_init_raf_all(fid, rafinfo, folder, rcts, cats):
	fid.write(folder + "\n")
	for i,r in enumerate(rafinfo[2]):
		if rcts[r,1] == 1: 
			str2w = str(int(rcts[r,2])) + " -> " + str(int(rcts[r,3])) + " + " + str(int(rcts[r,4])) + ", " + str(int(rafinfo[3][i])) + "\n"
			fid.write(str2w)
		else:
			str2w = str(int(rcts[r,3])) + " + " + str(int(rcts[r,4])) + " -> " + str(int(rcts[r,2])) + ", " + str(int(rafinfo[3][i])) + "\n"
			fid.write(str2w)
	fid.write("\n\n")

def write_acsCatalysis_file(path_file, catStr):
	print "\t|- save Catalysis file"
	fileName = os.path.join(path_file, '_acscatalysis.csv')
	np.savetxt(fileName, catStr, delimiter='\t', fmt='%d\t%d\t%d\t%d\t%f\t%f\t%f\t%d')

def write_acsReactions_file(path_file, rctStr, spontRatio=None, kspontass=None, kspontdiss=None):
	print "\t|- save Reaction file"
	if spontRatio and spontRatio > 0:
		for rct in rctStr:
			if ran.random() < spontRatio:
				if rct[1] == 1:rct[7] = kspontdiss
				else: rct[7] = kspontass
		
	fileName = os.path.join(path_file, '_acsreactions.csv')
	np.savetxt(fileName, rctStr, delimiter='\t', fmt='%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.10f')

def write_and_createInfluxFile(path_file, tmpFood):
	print "\t|- Save influx file"
	fname_influx = os.path.join(path_file, '_acsinflux.csv')
	fid_influx = open(fname_influx, 'w')
	prob = 1. / len(tmpFood)
	for id, f in enumerate(tmpFood):
		str2w = str(f) + "\t" + str(prob) + "\n"
		fid_influx.write(str2w)   
	fid_influx.close()

def write_and_create_std_nrgFile(path_file):
	print "\t|- Save energy file"
	fname_nrg = os.path.join(path_file, '_acsnrgbooleanfunctions.csv')
	fid_nrg = open(fname_nrg, 'w')
	str2w = "239\t1\n"
	fid_nrg.write(str2w)   
	fid_nrg.close()
	
def writeAllFilesAndCreateResFolder(pathFile, resFolderName, cats, rcts, food, spontRatio=None, kspontass=None, kspontdiss=None, conf=False):
	
	print "|- Save Structures"
	
	newdirAllResults = os.path.join(pathFile, resFolderName)
	if not os.path.isdir(newdirAllResults):
		try: os.mkdir(newdirAllResults)
		except: print "Impossible to create statistic directory", newdirAllResults; sys.exit(1)
		
	if conf == True: write_acsms_file(pathFile, conf)
	write_and_createInfluxFile(pathFile, food)
	write_and_create_std_nrgFile(pathFile)
	write_acsCatalysis_file(pathFile, cats)
	write_acsReactions_file(pathFile, rcts, spontRatio, kspontass, kspontdiss)
	
def write_acsms_file(path_file, nGen=10, nSim=1,nSec=1000,nRct=200000000,nH=0,nA=0,rs=0,dl=0,tssi=10,ftsi=0,nspmt=1,lfds=13,oc=0.0001,ecc=0,alf="AB",v=1e-18,\
					 vg=0,sd=0,nrg=0,rse=0,ncml=2,P=0.00103306,cp=0.5,mrevrct=0,rr=0,rrr=0,sr=0,K_ass=50,K_diss=25,K_cpx=50,K_cpxDiss=1,K_nrg=0,K_nrg_decay=0,K_spont_ass=0,\
					 K_spont_diss=0,moleculeDecay_KineticConstant=0.02,diffusion_contribute=0,solubility_threshold=0,influx_rate=0,maxLOut=3):
	print "\t|- save conf file"
	fname_conf = os.path.join(path_file, 'acsm2s.conf')
	fid_conf = open(fname_conf, 'w')
	
	fid_conf.write("# =================\n# ACSM2S PARAMETERS\n# =================\n\n# -----------------\n# SYSTEM PARAMETERS\n# -----------------\n\n")
	str2w = "nGEN="+str(nGen)+"\n"; fid_conf.write(str2w);
	str2w = "nSIM="+str(nSim)+"\n"; fid_conf.write(str2w);
	str2w = "nSeconds="+str(nSec)+"\n"; fid_conf.write(str2w);
	str2w = "nReactions="+str(nRct)+"\n"; fid_conf.write(str2w);
	str2w = "nHours="+str(nH)+"\n"; fid_conf.write(str2w);
	str2w = "nAttempts="+str(nA)+"\n"; fid_conf.write(str2w);
	str2w = "randomSeed="+str(rs)+"\n"; fid_conf.write(str2w);
	str2w = "debugLevel="+str(dl)+"\n"; fid_conf.write(str2w);
	str2w = "timeStructuresSavingInterval="+str(tssi)+"\n"; fid_conf.write(str2w);
	str2w = "fileTimesSaveInterval="+str(ftsi)+"\n"; fid_conf.write(str2w);
	str2w = "newSpeciesProbMinThreshold="+str(nspmt)+"\n"; fid_conf.write(str2w);
	str2w = "lastFiringDiskSpeciesID="+str(lfds)+"\n"; fid_conf.write(str2w);
	str2w = "overallConcentration="+str(oc)+"\n"; fid_conf.write(str2w);
	str2w = "ECConcentration="+str(ecc)+"\n"; fid_conf.write(str2w);
	str2w = "alphabet="+alf+"\n"; fid_conf.write(str2w);
	str2w = "volume="+str(v)+"\n"; fid_conf.write(str2w);
	str2w = "volumeGrowth="+str(vg)+"\n"; fid_conf.write(str2w);
	str2w = "stochDivision="+str(sd)+"\n"; fid_conf.write(str2w);
	str2w = "energy="+str(nrg)+"\n"; fid_conf.write(str2w);
	str2w = "ratioSpeciesEnergizable="+str(rse)+"\n"; fid_conf.write(str2w);
	str2w = "nonCatalyticMaxLength="+str(ncml)+"\n"; fid_conf.write(str2w);
	str2w = "reactionProbability="+str(P)+"\n"; fid_conf.write(str2w);
	str2w = "cleavageProbability="+str(cp)+"\n"; fid_conf.write(str2w);
	str2w = "main_rev_rct_allowed="+str(mrevrct)+"\n"; fid_conf.write(str2w);
	str2w = "reverseReactions="+str(rr)+"\n"; fid_conf.write(str2w);
	str2w = "revRctRatio="+str(rrr)+"\n"; fid_conf.write(str2w);
	str2w = "spontRct="+str(sr)+"\n"; fid_conf.write(str2w);
	str2w = "K_ass="+str(K_ass)+"\n"; fid_conf.write(str2w);
	str2w = "K_diss="+str(K_diss)+"\n"; fid_conf.write(str2w);
	str2w = "K_cpx="+str(K_cpx)+"\n"; fid_conf.write(str2w);
	str2w = "K_cpxDiss="+str(K_cpxDiss)+"\n"; fid_conf.write(str2w);
	str2w = "K_nrg="+str(K_nrg)+"\n"; fid_conf.write(str2w);
	str2w = "K_nrg_decay="+str(K_nrg_decay)+"\n"; fid_conf.write(str2w);
	str2w = "K_spont_ass="+str(K_spont_ass)+"\n"; fid_conf.write(str2w);
	str2w = "K_spont_diss="+str(K_spont_diss)+"\n"; fid_conf.write(str2w);
	str2w = "moleculeDecay_KineticConstant="+str(moleculeDecay_KineticConstant)+"\n"; fid_conf.write(str2w);
	str2w = "diffusion_contribute="+str(solubility_threshold)+"\n"; fid_conf.write(str2w);
	str2w = "solubility_threshold="+str(solubility_threshold)+"\n"; fid_conf.write(str2w);
	str2w = "influx_rate="+str(influx_rate)+"\n"; fid_conf.write(str2w);
	str2w = "maxLOut="+str(maxLOut)+"\n"; fid_conf.write(str2w);
 
	fid_conf.close()
