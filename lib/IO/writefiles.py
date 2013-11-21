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