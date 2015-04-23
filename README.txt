CARNESS APP :: Analysis Python Package
======================================

CARNESS APP is a python scripts package containing useful scripts to study and analyse the simulations carried out by means of the CARNESS simulator. 

System requirements
-------------------

Python Libraries: Numpy, Scipy, networkX, matplotlib 

Scripts
-------

- acsAttractorAnalysis.py -> Comparison between final states 
- acsAttractorAnalysisInTime.py -> Comparison between distinct states in time
- acsSCCanalysis.py -> Strongly connected components analysis
- acsSpeciesActivity.py -> Count how many times a species is a substrate, a catalyst or a product
- acsStateAnalysis.py -> Comparison between different times in the same run. Agglomerative statistics are computed as well
- acsDynStatInTime -> Trajectory in concentrations of each species are detected
- init.py system initializator
- prepareNewSim.py -> prepare the initial conditions of a new simulation starting from the final conditions of an other sim. 


controllare creazione di reazioni spontanee, 