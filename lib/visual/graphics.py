#!/usr/bin/python
# -*- coding: latin-1 -*-

'''
	This libraries contained the graphical functions useful for the CRN analysis
'''
import sys, os # Standard library
import numpy as np
import random as ran
from matplotlib import rc
from matplotlib import use
#use('Agg')
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from pylab import *
import networkx as nx
from networkx.algorithms import bipartite

# PLOT DEAFAULT VALUES
font = {'family' : 'serif',
        'weight' : 'normal',
        'size'   : 10,
        }
plt.rc('font', **font)
#plt.rc('text', usetex = True)
colr = ['k','r','b','g','c','m','y','burlywood','brown']

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

def PlotMatrix(tmpFilename, tmpX, tmpY, tmpXlabel, tmpYlabel, tmpLegend, tmpimgformat='png',tmpdpi=150):
	''' PlotMatrix is a function for plotting one or more lines.
		- tmpFilename is the path + file name of the file to save
		- tmpX can be a matrix column or a vector
		- tmpY can be a matrix or a vector
		- tmpTitle, tmpXlabel and tmpYlabel are a strings
		- tmpLegend can be a set of strings or can be set to None value if you dont't need legend 
	'''
	for i in range(tmpY.shape[1]):
		plt.plot(tmpX, tmpY[:,i], 'k',color=[ran.random(),ran.random(),ran.random()])
	plt.xlim(0,max(tmpX))
	plt.ylim(0,tmpY.max())
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
	filename = tmpFilename + '.' + tmpimgformat
	plt.savefig(filename, format=tmpimgformat, dpi=tmpdpi)
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

def plotBipartiteGraph(rcts, cats, tmpFolder, tmpFilename, imgname='_bipartite.pdf', savegraphimage=False, par_node_size=300, par_font_size=12):
	''' Function to create bipartite graphs starting from ctrs and cats structures

		:param rcts: reactions structures
		:param catalysis: catalysis structures
		:param tmpFolder: Save Folder	
	'''
	# Create Bipartite Complete Graph
	# Add catalyst -> reaction edges
	BIG = nx.DiGraph()
	# Add species from catalysis to graph
	BIG.add_nodes_from(cats[:,1], bipartite=0)
	# Add reactions to graph
	strrcts = [str(int(i)) for i in cats[:,2]]
	BIG.add_nodes_from(strrcts, bipartite=1) 
	# Add catalyst->reaction edges
	loe = [(x,strrcts[i]) for i,x in enumerate(cats[:,1])]
	BIG.add_edges_from(loe, weight=0.5)
	# Add reaction -> product nodes and edges
	# Add species from reactions to graph
	BIG.add_nodes_from(rcts[:,2], bipartite=0)
	BIG.add_nodes_from(rcts[:,3], bipartite=0)
	BIG.add_nodes_from(rcts[:,4], bipartite=0)
	# Create substrate reactions edges and reactions products edges
	newedges = []
	warningedges = []
	for r in rcts:
		if r[1]==0:
			 
			newedges.append((str(int(r[0])),r[2]))

			if nx.has_path(BIG,r[3],str(int(r[0]))): warningedges.append((r[3],str(int(r[0]))))
			else: newedges.append((r[3],str(int(r[0]))))
			
			if nx.has_path(BIG,r[4],str(int(r[0]))): warningedges.append((r[4],str(int(r[0]))))
			else: newedges.append((r[4],str(int(r[0]))))
		else:
			newedges.append((str(int(r[0])),r[3]))
			newedges.append((str(int(r[0])),r[4]))

			if nx.has_path(BIG,r[2],str(int(r[0]))): warningedges.append((r[2],str(int(r[0]))))
			else: newedges.append((r[2],str(int(r[0]))))
	
	BIG.add_edges_from(newedges, weight=1.0)
	BIG.add_edges_from(warningedges, weight=1.5)
	
	pos=nx.spring_layout(BIG)
	species_nodes, reactions_nodes = bipartite.sets(BIG)
	
	if savegraphimage:
		nx.draw_networkx_nodes(BIG,pos,nodelist=list(species_nodes), node_size=par_node_size, node_shape='o')
		nx.draw_networkx_nodes(BIG,pos,nodelist=list(reactions_nodes), node_size=par_node_size, node_shape='s',node_color="green")
		catalysis=[(u,v) for (u,v,d) in BIG.edges(data=True) if d['weight'] == 0.5]
		reactions=[(u,v) for (u,v,d) in BIG.edges(data=True) if d['weight'] == 1.0]
		warnings=[(u,v) for (u,v,d) in BIG.edges(data=True) if d['weight'] == 1.5]
		nx.draw_networkx_edges(BIG,pos,edgelist=reactions, width=0.5, alpha=0.5)
		nx.draw_networkx_edges(BIG,pos,edgelist=catalysis, width=0.5, alpha=0.5,edge_color='b',style='dashed')
		nx.draw_networkx_edges(BIG,pos,edgelist=warnings, width=1.5, edge_color='r')
		nx.draw_networkx_labels(BIG,pos, font_size=par_font_size)
		plt.axis('off')
		plt.savefig(os.path.join(tmpFolder, imgname))
		#plt.show()

	#nx.write_graphml(BIG, os.path.join(tmpFolder, tmpFilename))
	nx.write_pajek(BIG, os.path.join(tmpFolder, tmpFilename))
	plt.close()

def plotGraph(BIG, tmpFolder, tmpFilename, imgname='_bipartite.pdf', savegraphimage=False):
	'''
		Plot normal graphs
	'''

	if savegraphimage:
		pos=nx.random_layout(BIG)
		nx.draw_networkx_nodes(BIG,pos,node_shape='o')
		nx.draw_networkx_edges(BIG,pos, width=0.5)
		nx.draw_networkx_labels(BIG,pos)
		plt.axis('off')	
		#plt.show()
		plt.savefig(os.path.join(tmpFolder, imgname))

	nx.write_pajek(BIG, os.path.join(tmpFolder, tmpFilename))

