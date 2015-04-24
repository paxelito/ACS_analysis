CARNESS APP :: Analysis Python Package
======================================

CARNESS APP is a python scripts package containing useful scripts to study and analyse the simulations carried out by means of the CARNESS simulator. 

System requirements
-------------------

Python Libraries: python-numpy python-scipy python-matplotlib python-setuptools python-networkx python-deap graphviz libgraphviz-dev pkg-config python-pygraphvizNumpy, Scipy, networkX, matplotlib

All Scripts
-----------

* initializator.py 				-> system initializator
* graph_chemistry_analysis		-> artificial chemistry graph analysis
* acsDuplicationAnalysis		-> protcell duplication analysis
* acsAttractorAnalysis.py 		-> Comparison between final states 
* acsAttractorAnalysisInTime.py -> Comparison between distinct states in time
* acsSCCanalysis.py 			-> Strongly connected components analysis
* acsSpeciesActivity.py 		-> Count how many times a species is a substrate, a catalyst or a product
* acsStateAnalysis.py 			-> Comparison between different times in the same run. Agglomerative statistics are computed as well
* acsDynStatInTime 				-> Trajectory in concentrations of each species are detected
* prepareNewSim.py 				-> prepare the initial conditions of a new simulation starting from the final conditions of an other sim. 

Currently Documented
--------------------

* initializator.py 				-> system initializator
* graph_chemistry_analysis		-> artificial chemistry graph analysis
* acsDuplicationAnalysis		-> protcell duplication analysis

Release Notes
-------------

Version 001 - 20150424
**********************

* Documentation update 
* initializator.py :: the permission of the lunghers is not automatically set to 0755 



controllare creazione di reazioni spontanee, 