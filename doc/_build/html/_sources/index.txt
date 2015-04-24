.. CARNESS analysis python package

CARNESS analysis python package documentation
=============================================

The CARNESS analysis python package contains a set of python scripts devoted to 

1. the creation and the structural analysis of artificial catalytic reaction networks (CRS). 
2. the analysis of the dynamical outcomes of the artificial catalytic reaction networks simulated by means of the CAtalytic Reaction NEtworks Stochastic Simulator (`CARNESS  <https://github.com/paxelito/carness>`_).

System Requirements
-------------------

In order to use the scripts contained in this package, some further packages need to be installed. 
On UNIX systems you can install the necessary packages by the usual commands, for example on DEBIAN-based systems::

   $ sudo apt-get update
   $ sudo apt-get install python-numpy python-scipy python-matplotlib python-setuptools python-networkx python-deap graphviz libgraphviz-dev pkg-config python-pygraphviz

On WINDOWS system, I personally suggest that the simplest solution may be to rely on the alternative solutions proposed `here <https://www.python.org/download/alternatives/>`_.

Installing
----------
 
To get the packages (the package is entirely contained in a folder named "ACS_analysis") it is sufficent to clone it from GIT with the command::
   
   git clone https://github.com/paxelito/ACS_analysis.git

Since the package is still under development it is a good practice to check for possible updates or upgrades each time you use it. To seach for new updates please run the following command within the "ACS_analysis" folder ::

   git pull 

Usage
-----

To use the package you must run the opportune python script with the appropriate parameters. You can find a detailed description of the available scripts within this documentation. 
It is important to notice that not all the actually available scripts are depicted in this documentation. Hence, if you are familiar with python and with the CARNESS simulator, you can search into the package whether the desired analysis have been already implemented. 
Out objective is that of preparing the documentation for all the packages in the next few months. 


Current Available Scripts Analysis:

.. toctree::
   :maxdepth: 1

   initializator
   graph_chemistry_analysis
   acsDuplicationAnalysis
..   acsAttractorAnalysis
..  acsAttractorAnalysisInTime
..   acsBufferedFluxes
..   acsDynStatInTime
..   acsFromWim2Carness
..   acsSCCanalysis
..   acsSpeciesActivities
..   acsStatesAnalysis
..   lib
..   main
..   prepareNewSim
..   topology_analysis


Indices
=======

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

