#!/usr/bin/python
# -*- coding: latin-1 -*-

'''
	This libraries contained the graphical functions useful for the CRN analysis
'''
import sys, os # Standard library
import numpy as np
from matplotlib import rc
from matplotlib import use
use('Agg')
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from pylab import *
# PLOT DEAFAULT VALUES
font = {'family' : 'serif',
        'weight' : 'normal',
        'size'   : 10,
        }
plt.rc('font', **font)
#plt.rc('text', usetex = True)
colr = ['k','r','b','g','c','m','y']

#plt.rc('lines', linewidth=2)
#plt.rc("figure.subplot", left=(52/72.27)/width)
#plt.rc("figure.subplot", right=(width-10/72.27)/width)
#plt.rc("figure.subplot", bottom=(14/72.27)/height)
#plt.rc("figure.subplot", top=(height-7/72.27)/height)

def plot_np_txt_selcols_and_rows_allvarinrows(tmpfn,tmpdata,cols,vars,varcond,tmpxlab,tmpylab,tmplegend,tmpylim=None,tmpxlim=None,pstore='../',\
											  tmpimgformat='png',tmpdpi=150):
	''' Plot a numpy matrix with all the variables in rows. Variables are defined by the value of a particular colum
	'''
	tmpc = 0
	plt.figure(figsize=(4, 4))
	for i in vars:
		plt.plot(tmpdata[tmpdata[:,varcond]==i,cols[0]], tmpdata[tmpdata[:,varcond]==i,cols[1]], 'k', color=colr[tmpc])
		tmpc += 1
	
	plt.xlabel(tmpxlab, fontdict=font)
	plt.ylabel(tmpylab, fontdict=font)
	plt.legend(tmplegend, loc = 'best')
	plt.grid()
	if tmpxlim: plt.xlim(tmpxlim)
	if tmpylim: plt.ylim(tmpylim)
	
	actfolder = os.getcwd()
	os.chdir(pstore)	
	filename = tmpfn + '.' + tmpimgformat
	if os.path.isfile(filename): os.remove(filename) #remove the file before to create it (if exists)
	#plt.savefig(filename, format=tmpimgformat, dpi=tmpdpi)
	plt.savefig(filename, format=tmpimgformat)
	plt.close()
	os.chdir(actfolder)

def PlotMatrix(tmpFilename, tmpX, tmpY, tmpXlabel, tmpYlabel, tmpLegend):
	''' PlotMatrix is a function for plotting one or more lines.
		- tmpFilename is the path + file name of the file to save
		- tmpX can be a matrix column or a vector
		- tmpY can be a matrix or a vector
		- tmpTitle, tmpXlabel and tmpYlabel are a strings
		- tmpLegend can be a set of strings or can be set to None value if you dont't need legend 
	'''
	plt.plot(tmpX, tmpY, 'k')
	plt.xlim(0,max(tmpX))
	plt.ylim(0,max(tmpY))
#  	ax = plt.gca()
#  	for tick in ax.xaxis.get_major_ticks():
#  		tick.label1.set_fontsize(actfontsize)
#  		tick.label1.set_fontname('Times')
#  	for tick in ax.yaxis.get_major_ticks():
#  		tick.label1.set_fontsize(actfontsize)
#  		tick.label1.set_fontname('Times')

	plt.xlabel(tmpXlabel)
	plt.ylabel(tmpYlabel)
	plt.grid()
	if tmpLegend != None:
		plt.legend(tmpLegend, loc = 'best')
		
	plt.savefig(tmpFilename, format='eps', dpi=300)
	plt.close()
	#os.rename(tmpFilename, os.path.join('__10_stastisticFiles',tmpFilename))
	
def PlotHist(tmpFilename, tmpX, tmpXlabel, tmpYlabel,tmpimgformat='png',tmpdpi=150,num_bins=10):
	''' PlotMatrix is a function for plotting one or more lines.
		- tmpFilename is the path + file name of the file to save
		- tmpX can be a matrix column or a vector
	'''	
	# the histogram of the data
	n, bins, patches = plt.hist(tmpX, num_bins, normed=1, facecolor='green', alpha=0.5)
	# add a 'best fit' line
	# add a 'best fit' line
	y = mlab.normpdf(bins, np.mean(tmpX), np.std(tmpX))
	plt.plot(bins, y, 'r--')
	plt.xlabel(tmpXlabel)
	plt.ylabel(tmpYlabel)
	
	# Tweak spacing to prevent clipping of ylabel
	#plt.subplots_adjust(left=0.15)
	plt.savefig(tmpFilename, format=tmpimgformat, dpi=tmpdpi)
	plt.close()
	#os.rename(tmpFilename, os.path.join('__10_stastisticFiles',tmpFilename))