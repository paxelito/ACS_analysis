/** \mainpage Catalytic Rections Network Stochastic Simulator - CaRNeSS 4.8 (20131125.61)
 * \author Alessandro Filisetti
 * \version 4.8 (20131125.61)
 * \date 2013-11-25
 * sourceforge repository -- https://carness.svn.sourceforge.net/svnroot/carness/
 * git repository -- https://github.com/paxelito/carness
 *
 *
* This manual is divided in the following sections:
*- \subpage intro
*- \subpage pageInitStr
*- \subpage pageoutcomes
*- \subpage pageGillespie
*- \subpage pageInitializator
*
 * \page intro Essential informations
 * <hr>
 * The <b>Catalytic Reactions Network Stochastic Simulator (CaRNeSS)</b> is a computational model devoted to the simulation of theoretical complex catalytic networks composed of different
 * interacting molecular species. 
 * The model takes inspiration from the original model proposed by Stuart Kauffman in 1986, and describes systems composed of molecular species interacting by means of two possible reactions only, cleavage and condensation. One polymer is divided into two short polymers in the former case while two polymers are glued together forming a longer polymer in the latter case.
 * Each reaction must be catalyzed by another species in the system to occur, and one of the assumptions is that any chemical has an independent probability to catalyze a randomly chosen reaction.
 * Since the version 4.8 (20131026.60) spontaneous reactions are considered too.
 * It is important to notice that there are not indications about the chemical nature of the molecules,
 * species "A" may be both a polipeptide, an amminoacid, a particular protein domain or an RNA strenght.<br><br>
 *\section secUsage Using the simulator
 * To run the simulator open a terminal shell and type:<br><br>
 *	\c$path/executiveFile \c<configuration_File_Folder> \c<output_folder> \c<reaction_structures_folder><br>
 *		Examples:
 *      - Unix Based Systems:\c ~/Documents/project/acsm2s \c ~/Documents/.../confFileFolder/ \c ~/Documents/.../resFolder/ \c ~/Documents/.../StructuresFolder/
 *      - Win Systems: \c C:\\Documents\\project\\acsm2s.exe \c C:\\Documents\\...\\confFileFolder\\ \c C:\\Documents\\...\\resFolder\\ \c C:\\Documents\\...\\StructuresFolder\\
 *
 * <br><br>
 * \section sysreq System Requirement
 * <hr>
 * In order to have the simulator run correctly the recommended staff is reported:
 *		- MacOsX 10.4 or later, Linux (or in general a system UNIX based) or Windows OS (tests have been performed on Win7 and win Vista) as well
 *		- GCC (G++) compiler, or similar, installed (if you need to compile the software on your machine)
 *      - On MacOs system compile using \c g++ -Wall -ansi -lm -o CaRNeSS *.cpp
 *<br><br>
 * \section parameters Input Parameters :: acsm2s.conf
 * <hr>
 * All the system parameters are stored in a file called <b>acsm2s.conf</b>. Anyone can create his own configuration file paying attention to put "=" char between the parameter name and the the parameter value (NO SPACE BETWEEN THEM).<br>
 * Notice that the simulator does not create the initial structures but it simply loads the structures created by an external software and process them. Nevertheless the configuration file is fundamental to supply all the parameters to the simulation (during the simulation new entities may be created).
 * The simulator is provided with a structures initializator developed in MATLAB language by the group (a description of the initializator is provided in the main file "start.m") in which all the parameters we are going to describe are used to create the initial structures.
 * All parameters are reported below divided in three categories:
 *
 *		- System
 *		- Environment
 *		- Dynamic
 *
 *		Categories are useful only to help users in the parameter recognition within the configuration file. They are not handled from the software, if you like you can rearrange configuration file as you prefer, notice only that comments have to start with character <i>"#"</i>. Within the source code folder an example of the acsm2s.conf file is provided.<br>
 *      The following parameters are used both by the initializator and the simulator. Nvertheless it is ALWAYS necessary having a complete configuration file even if the structures have been already created.
 *		\subsection paramsystem System
 *		@param nGen (> 0) Number of generations. This parameter indicate how many times the simulation is stopped, concentration are set to the initial ones and the simulazion restart for other nSeconds seconds.
 *		@param nSIM (> 0) Number of simulations per generation starting with the same initial conditions (same data structures) but different random seed
 *		@param nSeconds (> 0) Number of seconds
 *		@param nReactions (> 0) Max number of reactions (the system will be stopped after nSeconds or after nReactions)
 *		@param randomSeed (>= 0) Random seed (if 0 the random seed is randomly created and the it is stored in the acsm2s.conf file saved in the results folder)
 *		@param nHours (>=0) Runtime limit (hours)
 *      @param nAttempts (>=0) Number of temptative in simulating the same network structure different random seed
 *      @param debugLevel (>= -1) Debug Level Runtime: different runTime message amounts (from -1 to 4, 0 is suggested)
 *		@param timeStructuresSavingInterval (> 0) All system structures (species, catalysis and reactions) are saved every <i>timeStructuresSavingInterval</i> seconds (simulation time)
 *      @param fileTimesSaveInterval (>= 0) Times data are stored in file times.csv every <i>fileTimesSaveInterval</i> seconds (If 0 reactions are stored continually)
 *		\subsection paramenv Environment
 *		@param newSpeciesProbMinThreshold (>=0) Minimal new species creation probability to allow system expansion
 *		@param lastFiringDiskSpeciesID (> 0) The ID of the last firing disk species.
 *		@param overallConcentration (> 0) The overall initial concentration that will be divided between all the initial species according to the selected initial distribution.
 *		@param ECConcentration (> 0) Incoming concentration of charged molecules per second.
 *		@param alphabet (string) Alphabet used in the simulation (e.g. <i>AGCT</i> for DNA, <i>ADEGFLYCWPHQIMTNKSRV</i> for proteins)
 *		@param volume (> 0) Volume of the container or protocell
 *		\subsection paramdyn Dynamic
 *		@param energy (0 or 1) 0 no energy in the system, 1 energy constraints are applied
 *      @param ratioSpeciesEnergizable (%) The probability for a species to be potentially energized by the energy carriers
 *		@param nonCatalyticMaxLength (>= 0) Max length of non catalytic species
 *		@param reactionProbability (from 0 to 1) Probability for a species to catalyze a reaction
 *		@param cleavageProbability (from 0 to 1) Cleavage probability (Condensation probability is 1 - cleavage probability)
 *		@param reverseReaction (0 or 1) Set to 1 to enable reverse reactions, 0 otherwise
 *      @param revRctRatio (>0) Ratio between forward and backward reactions, it is used in the creation of new reactions only (if reverseReactions = TRUE)
 *		@param spontRct (0 or 1) If 1 spontanoues reactions are considered, otherwise no
 *		@param K_ass (>= 0) Final Condensation kinetic constant (C.A + B --> AB + C) where A.C is the molecular complex composed of C (the catalyst) and A (the first substrate)
 *		@param K_diss (>= 0) Cleavage kinetic constant (AB --> A + B)
 *		@param K_cpx (>= 0) Complex formation kinetic constant (A + C(catalyst) --> C.A)
 *		@param K_cpxDiss (>= 0) Complex Dissociation kinetic constant (C.A --> A + C)
 *		@param K_nrg (>= 0) species phosphorilation kinetic constant
 *		@param K_nrg_decay (>= 0) de-energization kinetic constant
 *		@param moleculeDecay_KineticConstant (> 0) Molecule decay (efflux) kinetic Constant (Disregarded if the system is closed)
 *		@param K_spont_ass (>= 0) If spontaneous reactions are turned on this is the default kinetic rate for spontanoues condensations
 *		@param K_spont_diss (>= 0) If spontaneous reactions are turned on this is the default kinetic rate for spontanoues cleavages
 *		@param influx_rate (>= 0) Concentration per seconds (The species to insert in the system will be randomly chosen according to the _acsinflux.csv file). If equal to 0 the system is closed (maxLOut=0) or only the species that can cross the membrane come in and go out (maxLOut>0).
 *		@param maxLOut Maximum lenght of the species involved in the efflux process (\c influx_rate  > 0), equal to 0 indicates that all the species can be involved in the efflux process (no filter). If influx_rate = 0 the parameter indicates the species that can cross the semipermeable membrane of the protocell. <b>THE COUPLING BETWEEN INFLUX_RATE AND MAXLOUT INDICATES IF WE ARE SIMULATING A PROTOCELL OR A FLOW REACTOR</b>
 *		@param diffusion_contribute (KD) (0 or 0.5) if set to 0.5 the speed of molecules goes with the inverse of the square of the length, L^{-KD}
 *		@param solubility_threshold (> 0) Solubility Threshold, all the species longer than solubility_threshold precipitate
 *
 *<br><br>
 * \section Acknowledgments
 * <hr>
 *- University of Bologna, Interdepartment of industrial research (C.I.R.I)
 *- European Centre for Living Technology http://www.ecltech.org/
 *- Fondazione Venezia http://www.fondazionevenezia.it
 *- Alex Graudenzi to take care of the initilizator.
 *- Chiara Damiani to contribute to the development of the software.
 *- Roberto Serra, Marco Villani, Timoteo Carletti, Norman Packard, Ruedi Fuchslin and Stuart Kauffman for the essential hints.
 *- http://www.bedaux.net/mtrand/ for the pseudo-random Marseinne-Twister library for C++.
 *- http://perso.wanadoo.es/antlarr/otherapps.html for the poisson distribution generator numbers (acs_longInt random_poisson(acs_double tmpLambda, MTRand& tmpRandomGenerator)).
 *- Dr. Luca Ansaloni (luca.ansaloni@unimore.it) for the support but especially for the file handling functions and new Python development.
 *
 * \page pageInitStr Initial Data Structures
 * In order to proceed with the simulation the system needs from <b>4</b> to <b>6 (open system and energy on)</b> initial data structures files (an example for each file is located into the source code folder):
 *
 *      - \c acsm2s.conf (described in the \ref parameters section)
 *      - \c _acsspecies.csv - This file contains all the initial <b>species</b> with their proprieties
 *      - \c _acsreactions.csv - This file contains all the initial <b>reactions</b> with their proprieties
 *      - \c _acscatalysis.csv - This file contains all the <b>correspondances between species and reactions</b> with their proprieties
 *
 * and, if the system is open
 *
 *      - \c _acsinflux.csv - This file contains all the species belonging to the incoming flux
 *
 * and. if the system is energy based
 *
 *		- \c _acsnrgbooleanfunctions.csv - This file contains all the possible boolean functions associated with the reactions
 *
 * <br>
 *	\section subSpecies _acsspecies.csv
 * Columns description (each field is delimited using "\t"):
 *	  <table>
 *		<tr>
 *			<td>Identificator(ID)</td><td>Sequence</td><td>Concentration</td><td>Diffusion enhancement</td><td>Precipitation flag</td><td>Complex Dissociation Kinetic Constant</td><td>Complex Binding Point</td>
 *          <td>Already evaluated flag</td><td>Species Age</td><td>Number of reborns</td><td>Catalyst ID</td><td>Complex ID</td><td>Phosphorilation Kinetic constant</td><td>Charged Molecules Concentration</td><td>Concentration locked</td>
 *		</tr>
 *	  </table>
 *		- <i>ID</i>: Species Index
 *		- <i>Sequence</i>: Species sequence
 *		- <i>Concentration</i>: Total concentration of the species
 *		- <i>Diffusion enhancement</i>: Diffusion enhancement
 *		- <i>Precipitation flag</i>: If 0 species is precipited and each new molecules of this species will be precipited
 *		- <i>Dissociation Kinetic Constant</i>: Complex dissociation kinetic constant
 *		- <i>Binding point</i>: If the species is a complex this field indicates the division point between catalyst and substrate
  *		- <i>Evaluated</i>: If 1 the species is not virtual and all the reactions it catalyzes are created, if 0 the species is only potentially created
 *		- <i>Species Age</i>: Age (in seconds) of the species since its last cretion, each time that a species amount pass from 0 to > 1 the counter is resetted
 *		- <i>Number of reborns</i>: Number of times that a species amount pass from 0 to > 0
 *		- <i>Catalyst ID</i>: If the species is a complex (e.g. catalyst C forming a complex C.A with the substrate A) this is the CATALYST ID
 *		- <i>Substrate ID</i>: If the species is a complex (e.g. catalyst C forming a complex C.A with the substrate A) this is the SUBSTATE ID
 *		- <i>Phosphorilation Kinetic constant</i>: NOT USED NOW!!!
 *		- <i>Charged Molecules Concentration</i>: Concentration of the charged molecules belonging to the species.
 *		- <i>Concentration locked</i>: 0 -> Concentration of the species changes according to the reactions affecting it, 1 -> The concentration of the species is locked (permeable species, the concentration of the species is assumed to be constant within the cell and in the environment)
 *<br>
 *	\section subReactions _acsreactions.csv
 *	Columns description (each field is delimited using "\t"):
 *	  <table>
 *		<tr>
 *			<td>Identificator</td><td>Reaction type</td><td>Species 1</td><td>Species 2</td><td>Species 3</td>
 *			<td>Reaction counter</td><td>Energy type</td><td>k spont</td>
 *		</tr>
 *	  </table>
 *		- <i>Identificator</i>: Reaction ID
 *		- <i>Reaction type</i>: 0 Condensation, 1 Cleavage
 *		- <i>Species 1</i>: Product ID if reaction type = 1, Substrate ID reaction type = 0
 *		- <i>Species 2</i>: Product ID if reaction type = 0, Substrate ID reaction type = 1
 *		- <i>Species 3</i>: Product ID if reaction type = 0, Substrate ID reaction type = 1
 *		- <i>Reaction counter</i>:  Reaction occurrance counter
 *		- <i>Energy Type</i>: The reaction energetic configuration, 1 for endoergonic 0 for esoergonic
 *		- <i>Spontaneous Constant</i>: Spontaneous reaction constant
 *<br>
 *	\section subCatalysis _acscatalysis.csv
 * Columns description (each field is delimited using "\t"):
 *	  <table>
 *		<tr>
 *			<td>Identificator</td><td>Catalyst ID</td><td>Reaction ID</td><td>Catalysis counter</td><td>K condensation</td><td>K cleavage</td><td>K Complex Association</td><td>Complex creation substrate target</td>
 *		</tr>
 *	  </table>
 *		- <i>Identificator</i>: Catalysis ID
 *		- <i>Catalyst ID</i>: species (as catalyst) ID
 *		- <i>Reaction ID</i>: Reaction ID
 *		- <i>Catalysis counter</i>: Catalysis counter
 *		- <i>K final step end condensation kinetic constant</i>: Final step end condensation kinetic constant
 *		- <i>K Cleavage</i>: Cleavage Kinetic constant
 *		- <i>K complex association</i>: Complex association kinetic constant
 *		- <i>Complex creation substrate target</i>: Which substrate will be involved in complex creation, 1 or 2
 *<br>
 *	\section subInflux _acsinflux.csv
 * Columns description (each field is delimited using "\t"):
 *	  <table>
 *		<tr>
 *			<td>Identificator</td><td>Probabilitity</td>
 *		</tr>
 *	  </table>
 *		- <i>Identificator</i>: Species ID
 *		- <i>Probabilitity</i>: Probability to be selected when a species has to be inserted into the system
 *<br>
 *	\section subEnergy _acsnrgbooleanfunctions.csv
 * Columns description (each field is delimited using "\t"):
 *	  <table>
 *		<tr>
 *			<td>Energetic Boolean Function (decimal form)</td><td>Probabilitity</td>
 *		</tr>
 *	  </table>
 *		- <i>Energetic Boolean Function (decimal form)</i>: Energetic Boolean Function decimal form, e.g. 1001011010100 = 4820
 *		- <i>Probabilitity</i>: Probability to be selected when a new reaction has to be created
 *
 * \page pageoutcomes Outcomes
 * During the simulation <b>8</b> different files concerning species, reactions and catalysis tables, times, reactions parameters, living species, living species amounts and concentration are produced. In the following you find
 * a brief description of all files and some indications about the different time storage.
 *	\section outspecies Species
 *	  Every environment#timeStructuresSavingInterval seconds species structures will be saved in files called species_<generationNumber>_<simulationNumber>_<reactionNumber>.csv
 *	  (e.g. <b>species_1_06_0016933.csv</b> means species structure, generation number 1, simulation number 6, reaction number 16933).<br>Columns description (each field is delimited using "\t") is the same described above, section \ref subSpecies .
 *
 *	\section outreactions Reactions
 *	Every environment#timeStructuresSavingInterval seconds reactions structures will be saved in files called reactions_<generationNumber>_<simulationNumber>_<reactionNumber>.csv
 *	  (e.g. <b>reactions_1_06_0016933.csv</b> means reactions structure, generation number 1, simulation number 6, reaction number 16933). <br>Columns description (each field is delimited using "\t") is the same described above, section \ref subReactions .
 *
 *	\section outcatalysis Catalysis
 *  Every environment#timeStructuresSavingInterval seconds catalysis structures will be saved in files called catalysis_<generationNumber>_<simulationNumber>_<reactionNumber>.csv
 *	  (e.g. <b>catalysis_1_06_0016933.csv</b> means catalysis structure, generation number 1, simulation number 6, reaction number 16933). <br>Columns description (each field is delimited using "\t")is the same described above, section \ref subCatalysis :
 *
 *	\section outTimes Times
 *	Every simulation generates a times file called times_<generationNumber>_<simulationNumber>.csv (e.g. <b>times_1_03.csv</b> means times file, generation 1, simulation 3) containing values concerning <b>each events</b>
 *  occurred during the simulation. This file is <b>created at the beginning of the simulation</b> and it is <b>updated after each event</b>, hence this file will contain a number of rows equal to the number of reaction (of any type) occurred
 *  during the simulation. <br>Columns description (each field is delimited using "\t"):
 *	  <table>
 *		<tr>
 *			<td>Identificator (reaction)</td><td>Reaction Time</td><td>Gillespie ID selected</td><td>Reaction Type</td><td>Number of possible reactions</td><td>Computational Time (sec)</td>
 *			<td>Number of Species</td><td>Number of Molecules</td><td>Number of complex species</td><td>Number of complexes</td>
 *                      <td>Number of bricks</td><td>Gillespie Computational Time</td><td>Reaction Process Computational Time</td><td>Various Processes Computational Time</td><td>New species creation probability</td>
 *		</tr>
 *	  </table>
 *		- <i>Identificator (reaction)</i>: Number of the reaction (STEP)
 *		- <i>Reaction Time</i>: Time of the reaction
 *		- <i>Gillespie ID selected</i>: Gillespie structure reaction ID
 *		- <i>Reaction Type</i>: 0 condensation, 1 cleavage, 2 complex formation, 3 complex degradation, 4 species efflux, 6 endergonic cleavage, 7 endoergonic condensation, 8 endoergonic complex creation
 *		- <i>Number of possible reactions</i>: Total number of possible reactions according to the Gillespie algorithm computation
 *		- <i>Computational Time (ms)</i>: Computational time between two successive reactions
 *		- <i>Number of Species</i>: Number of species with at least one molecule
 *		- <i>Number of Molecules</i>: Number of molecules
 *		- <i>Number of complex species</i>: Number of complex species with at least one molecule
 *		- <i>Number of complexes</i>: Number of complex molecules
 *		- <i>Number of bricks</i>: Number of single bricks (monomers) present in the system
 *		- <i>Gillespie Computational Time</i>: Number of milliseconds necessary to complete the Gillespie task
 *		- <i>Reaction Process Computational Time</i>: Number of milliseconds necessary to perform a reaction (and evaluate products) once that the Gillespie algorithm has selected the reaction
 *		- <i>Various Processes Computational Time</i>: Number of milliseconds necessary to perform several tasks not correlated with the simulation of the phenomena
 *		- <i>New species creation probability</i>: Given the state of the system, probability to create a new species
 *		- <i>Reverse Reaction Probability</i>: Given the state of the system, probability for a reverse reaction to occur
*
 *	\section outReactions Reactions_parameters
 *  Every simulation generates a reactions parameters file called reactions_parameters_<generationNumber>_<simulationNumber>.csv (e.g. <b>reactions_parameters_1_07.csv</b> means reactions parameters file, generation 1 simulation
 *  number 7) containing informations on <b>cleavage and condensation</b> reactions. This file is <b>generated at the beginning of each simulation</b> and it is <b>updated each time that a cleavage or a condensation occur</b>. 
 *  <br>Columns description (each field is delimited using "\t"):
 *	  <table>
 *		<tr>
 *			<td>Reaction ID</td><td>Reaction Time</td><td>Reaction Type</td><td>Catalyst ID</td><td>Species 1 ID</td><td>Species 2 ID</td><td>Species 3 ID</td>
 *			<td>Charged Molecules</td><td>Charged Concentration</td><td>Gillespie Score Average</td><td>Gillespie Score Standard Deviation</td><td>Entropy</td><td>New species creation probability</td>
 *		</tr>
 *	  </table>
 *		- <i>Reaction ID</i>: This ID referes to the reaction number (step)
 *		- <i>Reaction Time</i>: Time of the reaction
 *		- <i>Reaction Type ID</i>: 0 condensation, 1 cleavage, 6 endoergonic cleavage, 7 endoergonic condensation
 *		- <i>Catalyst ID</i>: Catalyst ID into the species table
 *		- <i>Species 1 ID</i>: product ID in case of condensation, substrate in case of cleavage
 *		- <i>Species 2 ID</i>: substrate ID in case of condensation, product in case of cleavage
 *		- <i>Species 3 ID</i>: substrate ID in case of condensation, product in case of cleavage
 *		- <i>Charged Molecules</i>: Number of charged molecules
 *		- <i>Charged Concentration</i>: Concentration of the charged molecules
 *		- <i>Gillespie Score Average</i>: Average of the all gillespie algorithm scores 
 *		- <i>Gillespie Score Standard Deviation</i>: Standard Deviation of the all gillespie algorithm scores
 *		- <i>Entropy</i>: Entropy Misure (test) based on the probabilities that reactions occur
 *		- <i>New species creation probability</i>: Given the state of the system, probability to create a new species
 *
 *
 *	\section outLivSpecies Living_species
 *  Each simulation generates a living species file called livingSpecies_<generationNumber>_<simulationNumber>.csv (e.g. <b>livingSpecies_1_07.csv</b> means living species file, generation 1 simulation 
 *  number 7) containing all the living species IDs (living species are those species with concentration greater than 0) for each reaction in each raw. This file is <b>generated at the beginning of each simulation</b> and it is <b>updated each time that a cleavage or a condensation occur</b>. 
 *  This file has a variable number of columns equal to the number of living species plus 2 (Reaction ID and time of the reaction)
 *  <br>Columns description (each field is delimited using "\t"):
 *	  <table>
 *		<tr>
 *			<td>Reaction ID</td><td>Reaction Time</td><td>living species IDs..., one for column</td>
 *		</tr>
 *	  </table>
 *		- <i>Reaction ID</i>: This ID referes to the reaction number (step)
 *		- <i>Reaction Time</i>: Time of the reaction
 *		- <i>Living Species IDs..., one for column</i>: Each column from the thirth to the last one contains living species ID
 *
 *	\section outLivAmounts Living_species_amounts
 *  Each simulation generates a living species amount file called livingAmount_<generationNumber>_<simulationNumber>.csv (e.g. <b>livingAmount_1_07.csv</b> means living species amount file, generation 1, simulation 
 *  number 7) containing all the living species amount (living species are those species with concentration greater than 0) for each reaction in each raw. <b>This file is tightly correlated with the <i>Living_species</i> file</b>. Each row contains the total amount of molecules belonging to the species
 *  indicated by its ID in the livingSpecies file. This file is <b>created at the beginning of each simulation</b> and it is <b>updated each time that a cleavage or a condensation occur</b>. 
 *  This file has a variable number of columns equals to the number of living species at time t
 *  <br>Columns description (each field is delimited using "\t"):
 *	  <table>
 *		<tr>
 *			<td>living species amounts..., one for column</td>
 *		</tr>
 *	  </table>
 *		- <i>Living Species amounts..., one for column</i>: Each values is related to the species indicated in the livingSpecies file. Reaction ID and time is that of the livingSpecies file
 *
 *	\section outLivConcentrations Living_species_concentration
 *  Each simulation generates a living species amount file called livingConcentration_<generationNumber>_<simulationNumber>.csv (e.g. <b>livingConcentration_1_07.csv</b> means living species amount file, generation 1, simulation 
 *  number 7) containing all the living species amount (living species are those species with concentration greater than 0) for each reaction in each raw. <b>This file is tightly correlated with the <i>Living_species</i> file</b>. Each row contains the total amount of molecules belonging to the species
 *  indicated by its ID in the livingSpecies file. This file is <b>created at the beginning of each simulation</b> and it is <b>updated each time that a cleavage or a condensation occur</b>. 
 *  This file has a variable number of columns equals to the number of living species at time t
 *  <br>Columns description (each field is delimited using "\t"):
 *	  <table>
 *		<tr>
 *			<td>living species amounts..., one for column</td>
 *		</tr>
 *	  </table>
 *		- <i>Living Species amounts..., one for column</i>: Each values is related to the species indicated in the livingSpecies file. Reaction ID and time is that of the livingSpecies file
 *
 * \page pageGillespie Gillespie Class
 * The Gillespie class contains the information concerning the Gillespie algorithm propensity with other useful information.
 * In particular it has 4 properties, namely mol_I, mol_II, mol_III, mol_IV containing the species ID's related to different species type according to the different admitted reactions.
 * In the following a table summarizing the different meanings for the different reactions is reported.
 * <table>
 * <tr>
 * 	<td>Reaction Type</td><td>Reaction Type Code</td><td>MOL 1</td><td>MOL 2</td><td>MOL 3</td><td>MOL 4</td>
 * </tr>
 *  <tr>
 * 	<td>Complex Formation</td><td>2</td><td>Catalyst</td><td>First Substrate</td><td>Catalysis ID</td><td>Second Substrate</td>
 * </tr>
 * <tr>
 * 	<td>Final Condensation</td><td>0</td><td>Catalyst</td><td>Substrate</td><td>Product</td><td>Complex</td>
 * </tr>
 * <tr>
 * 	<td>Cleavage</td><td>1</td><td>Substrate</td><td>Product 1</td><td>Product 2</td><td>Catalyst</td>
 * </tr>
 * <tr>
 * 	<td>Spontaneous Cleavage</td><td>11</td><td>Substrate</td><td>Product 1</td><td>Product 2</td><td>//</td>
 * </tr>
  * <tr>
 * 	<td>Spontaneous Condensation</td><td>10</td><td>Product</td><td>Substrate 1</td><td>Substrate 2</td><td>//</td>
 * </tr>
 * </table>
 *
 *
 * \page pageInitializator The initializator (a very brief description)
 * The initializator provided with the simulator is located in the \c initializator folder (within the source code folder) and it is developed in Matlab code. All the parameters are set in the \c start.m file (from line 22 to line 63). In addition you find the parameters related to the name of the folder that will contain the simulation (\c simFolder.name), the path where that folder will be created (\c simFolder.path) and the number of different network ensambles to create (\c simFolder.nets). It is important to notice that the initializator has not been thought to be shared, so it is not too much user friendly to be manipulated. Nevertheless it could be very useful with a little bit of practice.
 *  \section initSensitivity Screening Parameter
 * To initialize structures to perform a sensitivity analysis of a specific parameters follow the following instruction:
 *  - Change the name of the array at row 13 with the name of the parameter you want analyze. This array contains the values of the paramter. Remember to insert the old parameter you are changing in the parameters list with its single value, otherwise it would be a missing parameter.
 *  - In the array at row 14 (\c nome_folder) you have put the same number of elements of the array containing the screening values. This array contains a numeric tag of the values (used to create the simulations folder names) contained in the values array.
 *  - In the code of the \c start.m file (rows 122 to 150) you must change:
 *      - row 121: \c [rows,b]=size(reactionProbability); --> [rows,b]=size(name_of_the_parameter_you_have_trasformed_in_array);
 *      - row 146: \c inizializzatore_ACS(... reactionProbability(i), ...); -> inizializzatore_ACS(... name_of_the_parameter_you_have_trasformed_in_array(i), ...); (remember to remove the index to the reactionProbability parameter
 */


#include "acs_headers.h"
#include "environment.h"

void saveToFile(string tmpSavingPath, environment *tmpEnvironment, acs_int tmpGen, acs_int tmpSim, acs_int tmpStep);
void saveTimesToFile(string tmpSavingPath, environment *tmpEnvironment, acs_int tmpGen, acs_int tmpSim, acs_int tmpStep);
void saveInitialConditionsToFile(string tmpSavingPath, environment *tmpEnvironment, acs_int tmpGen, acs_int tmpSim, acs_int tmpStep);

int main (int argc, char *argv[]) {
	
    setlocale(LC_NUMERIC,"en_US"); //TR Check whether or not it can be removed when C++ standard
	
	MTRand rndDoubleGen; /// double random number generator
	
	/*---------------\
    |                |
	| INITIALIZATION |
    |       		 |
	\ --------------*/ 
	
		/*----------------
		CREATE ENVIRONMENT
		-----------------*/ 

        // Declare and inizialize the environment
        environment* puddle = new environment(argv[1]);
	
		//INITIZIALIZE PSEUDO RANDOM NUMBER GENERATOR
		rndDoubleGen.seed(puddle->getRandomSeed());
		
		/*----------------------------------
		CREATE INITIAL MOLECULE POPULATION
		----------------------------------*/
	
            //LOAD SPECIES FROM FILE
            if(!puddle->createInitialMoleculesPopulationFromFileSTD(argv[3]))
                ExitWithError("createInitialMoleculesPopulationFromFile", "Problem with the species STANDARD loading process");
            // LOAD INFLUX LAYERS FROM FILE (if the system is open with a simulated flux)
			if(puddle->getInflux() > 0)
			{
                if(!puddle->createInfluxLayersFromFileSTD(argv[3]))
					ExitWithError("CreateInfluxLayersFromFile", "Problem with influx layers loading process");
			}
                        // LOAD BOOLEAN FUNCTION CONCERNING THE ENERGY CONFIGURATION
			if(puddle->getEnergy() <= 1)
			{
                //if(!puddle->createNrgBooleanFunctionsFromFile(a.arguments().at(3)))
                if(!puddle->createNrgBooleanFunctionsFromFileSTD(argv[3]))
						ExitWithError("createNrgBooleanFunctionsFromFile", "Problem with rct Bool fncs loading process");
			}
           //LOAD REACTIONS STRUCTURE FROM FILE (standard C++ libraries)
           if(!puddle->createInitialReactionsLayerFromFileSTD(argv[3]))
				ExitWithError("createInitialReactionLayerFromFile", "Problem with the reactions loading process");

           //LOAD CATALYSIS STRUCTURE FROM FILE (catalysis links species with reactions catalyzed)
            if(!puddle->createInitialCatalysisLayerFromFileSTD(argv[3]))
				ExitWithError("createInitialCatalysisLayerFromFile", "Problem with the catalysis loading process");
		
		if(puddle->getDebugLevel() >= MEDIUM_DEBUG){puddle->printInitialCondition();}
	
                // SAVE TO FILE INITIAL CONDITIONS
        saveInitialConditionsToFile(argv[2], puddle, 0, 0, 0);
		
		if(!puddle->structureCoherenceCheckUp())
		{
			cout << endl;
            ExitWithError("structureCoherenceCheckUp", "PROBLEM WITH STRUCTURE COHERENCE... THE SIMULATION WILL BE ABORTED!!!");
		}
	
	
	// - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 
	
		/*------------\
		 |            |
		 | SIMULATION |
		 |            |
		 \ ----------*/
		
		//STORE INITIAL STRUCTURES DATA IN ORDER TO REINITIALIZE THE STRUCTURE AFTER EACH SIMULATION
		puddle->storeInitialStructures();
		
        clock_t tStart = clock();
		
		if(puddle->getDebugLevel() >= RUNNING_VERSION)
		{
			cout << endl;
			cout << "*** SIMULATION START ***" << endl;
			cout << endl;
		}
				
        // +-+-+-+-+-+-+-+-+-+-+-+-+-+-+
        // SIMULATIONS
        // +-+-+-+-+-+-+-+-+-+-+-+-+-+-+
		
        acs_int totalNumberOfSimulations = puddle->getNsim(); // Simulation counter
        acs_int totalNumberOfGenerations = puddle->getNgen(); // generation are not used now (it will be usefull with protocell division)
        acs_int actSTEP; // declare simulation step
        acs_int previousStepLastStructuresSaving; // declare previousStep
        acs_int previousStepLastTimesSaving; // declare previousStep

        for(acs_int actSIM = 1; actSIM <= puddle->getNsim(); actSIM++)
        {
            if(puddle->getDebugLevel() >= RUNNING_VERSION){cout << "\n|- - - - - - - - - - - - - - - " << actSIM << endl;}
            if(puddle->getDebugLevel() >= RUNNING_VERSION){cout << "|- SIMULATION NUMBER " << actSIM << endl;}
            if(puddle->getDebugLevel() >= RUNNING_VERSION){cout << "|- - - - - - - - - - - - - - - \n" << actSIM << endl;}

		
			// +-+-+-+-+-+-+-+-+-+-+-+-+-+-+
            // GENERATIONS PER SIMULATIONS
            // +-+-+-+-+-+-+-+-+-+-+-+-+-+-


            for(acs_int actGEN = 1; actGEN <= totalNumberOfGenerations; actGEN++)
            {
                // IF THE ATTEMPTS ARE MORE THAN THE MAX NUMBER THIS SIMULATION IS STOPPED
                if((puddle->getCurrentAttempts() <= puddle->getMAXattempts()) || (puddle->getMAXattempts() == 0))
                {
                    if(actGEN > 1)
                    {
                        // Reset structures stats and set concentrations to initial values
                        puddle->resetConcentrationToInitialConditions();
                        saveToFile(argv[2], puddle, actGEN, actSIM, 0);

                    }
                        acs_double dataStoredCounter = 0;
                        acs_double TimesStoredCounter = 0;
                        if(puddle->getDebugLevel() >= MINIMAL_PROMPT)
                        {
                             cout << "|- GEN NUMBER " << actGEN << " OF " << puddle->getNgen()
                             << " - SIM NUMBER " << actSIM << " OF " << totalNumberOfSimulations << endl;
                        }

                        tStart = clock();
                        actSTEP = 1;
                        previousStepLastStructuresSaving = 1;
                        previousStepLastTimesSaving = 1;

                        // +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
                        // SECONDS / REACTIONS PER SIMULATION
                        // +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

                        while((puddle->getActualTime() <= puddle->getNseconds()) & (actSTEP <= puddle->getNreactions())	)
                        {
                            // IF NUMBER OF MILLISECONDS IS LESS THAN THE MAX NUMBER
                            if(( (((float)clock() - tStart) / CLOCKS_PER_SEC) < (puddle->getMAXhours()*60*60)) || (puddle->getMAXhours() == 0))
                            {
                                //GILLESPIE COMPUTATION (CORE OF THE SOFTWARE)
                            	if(puddle->getDebugLevel() == SMALL_DEBUG)cout << "Step " << actSTEP << endl;
                            	try{
                                	if(!puddle->performOPTGillespieComputation(rndDoubleGen, tStart, actGEN, actSIM, actSTEP, argv[2]))
                                        	ExitWithError("performGillespieComputation", "Problems with the Gillespie computation");
                            	}catch(exception&e){
                           	     cout << "Source Code Line: " << __LINE__ << endl;
                           	     cerr << "exceptioncaught:" << e.what() << endl;
                           	     ExitWithError("MAIN function __ Gillespie computation","exceptionerrorthrown");
                            	}
                                if(puddle->getDebugLevel() == SMALL_DEBUG)
                                	puddle->showGillEngagementInSpecies();

                                // DISPLAY SIMULATION CONTROL VARIABLES
                                if(puddle->getDebugLevel() >= RUNNING_VERSION)
                                {
                                        if((actSTEP % 1000 == 0) || (actSTEP == 1) ||
                                           (puddle->getNseconds() - puddle->getActualTime() < 0.0001) ||
                                           (puddle->getMols() == puddle->getOverallLoadedMolsCounter()))
                                        {
                                                cout<< "--------------------------------------------------------------------" << endl
                                                    << "G: "<< actGEN << "/" << totalNumberOfGenerations
                                                    << " - S: " << actSIM << "/" << totalNumberOfSimulations
                                                    << " - T: " << puddle->getActualTime() << "/" << puddle->getNseconds()
                                                    << " - R: " << actSTEP // << "/" << puddle->getNreactions()
                                                    << " - CT (seconds): " << ((float)clock() - tStart) / CLOCKS_PER_SEC
                                                    << " - Gill: " << puddle->getNumberOfGillespieCOPYpossibleRcts() << endl
                                                    << "\t- ENVIRONMENT" << endl
                                                    << "\t\t|- S: " << puddle->getNspecies()
                                                    << " - NS: " << puddle->getNnewSpecies()
                                                    << " - M: " << puddle->getMols()
                                                    << " - NM: " << puddle->getNewMols() << endl
                                                    << "\t\t|- Cp: " << puddle->getNcpx()
                                                    << " - #Cp: " << puddle->getNcpxMols()
                                                    << " - Volume: " << puddle->getVolume() << endl
                                                    << "\t\t|- M: " << puddle->getTotalNumberOfMonomers()
                                                    << " - Mols+Complex: " << puddle->getMols()+puddle->getNcpxMols()
                                                    << " - Tot Conc: " << double(puddle->getMols()+puddle->getNcpxMols())/(puddle->getVolume()*AVO) << endl
                                                    << "\t- ENERGY" << endl
                                                    << "\t\t|- Loaded Mols: " << puddle->getOverallLoadedMolsCounter()
                                                    << " - Loaded (second Method): " << puddle->getTotNumberOfChargedMols() << endl
                                                    << "\t- GILLESPIE INFO" << endl
                                                    << "\t\t|- Gil Mean: " << puddle->getGillespieMean()
                                                    << " - Gil SD: " << puddle->getgillespieSD()
                                                    << " - Gil Entropy: " << puddle->getgillespieEntropy() << endl
                                                    << "\t\t|- Gil NS ratio: " << puddle->getRatioBetweenNewGillTotGill()
                                                    << " - Back and Forw Ratio: " << puddle->getRatioBetweenBackandForw() << endl
                                                    << "\t- REACTIONS COUNTERS" << endl
                                                    << "\t\t|- Theoretical Average Connectivity:" << (double)puddle->getNumberOfCatalysis() / puddle->getTotalNumberOfPossibleCatalysts() << endl
                                                    << "\t\t|- Cond: " << puddle->getCondensationCounter()
                                                    << " - Endo Cond: " << puddle->getEndoCondensationCounter() << endl
                                                    << "\t\t|- Cleav: " << puddle->getCleavageCounter()
                                                    << " - Endo Cleav: " << puddle->getEndoCleavageCounter() << endl
                                                    << "\t\t|- Spont Cleav: " << puddle->getSpontDissCounter()
                                                    << " - Spont Cond: "<< puddle->getSpontAssCounter() << endl
                                                	<< "\t\t|- Cpx: " << puddle->getCpxFormCounter()
                                                	<< " - Cpx Diss: " << puddle->getCpxDissCounter() << endl;
                                        }
                                }


                                // SAVE STRUCTURES TO FILE EVERY puddle->getTimeStructuresSavingInterval() SECONDS
                                if(puddle->getActualTime() >= (puddle->getTimeStructuresSavingInterval() + dataStoredCounter))
                                {
                                        saveToFile(argv[2], puddle, actGEN, actSIM, actSTEP);
                                        dataStoredCounter = dataStoredCounter + puddle->getTimeStructuresSavingInterval();                                        
                                        previousStepLastStructuresSaving = actSTEP; // save last step moment
                                }

                                // STORE ON FILE TIMES EVERY fileTimesSaveInterval seconds (if at least something happen)
                                if((puddle->getActualTime() >= (puddle->getFileTimesSavingInterval() + TimesStoredCounter)) ||
                                        (puddle->getActualTime() == 0) || (puddle->getFileTimesSavingInterval() == 0))
                                {
                                    //cout << puddle->getActualTime() << " " << puddle->getFileTimesSavingInterval() << " " << TimesStoredCounter << endl;
                                    saveTimesToFile(argv[2], puddle, actGEN, actSIM, actSTEP);
                                    if(puddle->getActualTime() > 0)
                                    {
                                        TimesStoredCounter = TimesStoredCounter + puddle->getFileTimesSavingInterval();
                                        previousStepLastTimesSaving = actSTEP;
                                    }
                                }

                                actSTEP++; // step update

                            }else{ // if(( (((float)clock() - tStart) / CLOCKS_PER_SEC) < (puddle->getMAXhours()*60*60)) || (puddle->getMAXhours() == 0))
                                // Increase number of attempts
                                puddle->increaseAttempts();
                                puddle->resetConcentrationToInitialConditions();
                                // STOP WHILE
                                actGEN = actGEN - 1;
                                break;

                            }
                        } // while((puddle->getActualTime() <= puddle->getNseconds()) & (actSTEP <= puddle->getNreactions()))

                        // Final species ages update
                        puddle->updateSpeciesAges();

                        // SAVE FINAL STRUCTURES TO FILE if actSTEP is different from the previous one
                        if(previousStepLastStructuresSaving < actSTEP-1)
                            saveToFile(argv[2], puddle, actGEN, actSIM, actSTEP);

                        if(previousStepLastTimesSaving < actSTEP-1)
                            saveTimesToFile(argv[2], puddle, actGEN, actSIM, actSTEP);

                        //CHECK STRUCTURES, THE COHERENCE OF THE INTERNAL STRUCTURES ARE CONTROLLED
                        if(!puddle->structureCoherenceCheckUp())
                        {
                            cout << endl;
                            ExitWithError("structureCoherenceCheckUp", "PROBLEM WITH STRCTURE COHERENCE... BE CARE TO THE SIMULATION OUTCOMES!!!");
                        }
                }else{ // for if(puddle->getCurrentAttempts() <= puddle->getMAXattempts())
                  // SIMULATION IS STOPPED
                  ExitWithError("MAIN", "FATAL ERROR: The number of MAXIMUM temptatives has been exceeded!!!");
                  break;

                }
            } // end for(acs_int actGEN = 1; actGEN <= puddle->getNgen(); actGEN++)
            // RESET STRUCTURES and TIME TO THE INITIAL VALUES FOR THE NEXT SIMULATION
            puddle->clearAllStructures();
        } // end for(acs_int actSIM = 1; actSIM <= puddle->getNsim(); actSIM++)
					
	//DELETE ALL MAIN HEAP OBJECTS
	delete puddle;
    return 0;
 }


/*---------------------------	
 MAIN SAVING FUNCTIONS
 ---------------------------*/

/**
 Save to file structures at step tmpStep
 @version 2.0
 @param string tmpSavingPath Saving files path
 @param environment *tmpEnvironment environment instance reference
 @param tmpSim Current simulation
 @param acs_int Current step
 @date 2013/07/13
 */
void saveToFile(string tmpSavingPath, environment *tmpEnvironment, acs_int tmpGen, acs_int tmpSim, acs_int tmpStep)
{
	// SAVE TO FILE INITIAL CONDITIONS
    if(!tmpEnvironment->saveSpeciesStructureSTD(tmpGen, tmpSim, tmpStep, tmpSavingPath))
		ExitWithError("saveSpeciesStructure", "PROBLEM WITH SAVING SPECIES STRUCTURE TO FILE");
    if(!tmpEnvironment->saveReactionsStructureSTD(tmpGen, tmpSim, tmpStep, tmpSavingPath))
		ExitWithError("saveReactionsStructure", "PROBLEM WITH SAVING REACTIONS STRUCTURE TO FILE");
    if(!tmpEnvironment->saveCatalysisStructureSTD(tmpGen, tmpSim, tmpStep, tmpSavingPath))
		ExitWithError("saveCatalysisStructure", "PROBLEM WITH SAVING CATALYSIS STRUCTURE TO FILE");	
}

/**
 Save TIMES to file
 @version 2.0
 @param string tmpSavingPath Saving files path
 @param environment *tmpEnvironment environment instance reference
 @param tmpSim Current simulation
 @param acs_int Current step
 @date 2013/07/03
 */
void saveTimesToFile(string tmpSavingPath, environment *tmpEnvironment, acs_int tmpGen, acs_int tmpSim, acs_int tmpStep)
{
    // SAVE TO FILE INITIAL CONDITIONS
    if(!tmpEnvironment->saveTimesSTD(tmpGen, tmpSim, tmpStep, tmpSavingPath))
        ExitWithError("saveTimes", "PROBLEM WITH SAVING TIMES");
}

/**
 Save to file all the INITIAL structures
 @version 2.0
 @param string tmpSavingPath Saving files path
 @param environment *tmpEnvironment environment instance reference
 @param tmpSim Current simulation
 @param acs_int Current step
 @date 2013/07/03
 */
void saveInitialConditionsToFile(string tmpSavingPath, environment *tmpEnvironment, acs_int tmpGen, acs_int tmpSim, acs_int tmpStep)
{
	// SAVE TO FILE INITIAL CONDITIONS
    if(!tmpEnvironment->saveConfigurationFileSTD(tmpSavingPath))
		ExitWithError("saveConfigurationFile", "PROBLEM WITH SAVING CONFIGURATION FILE");	
	if (tmpEnvironment->getInflux() > 0)
	{
        if(!tmpEnvironment->saveInfluxStructureSTD(tmpSavingPath))
			ExitWithError("saveInfluxStructure", "PROBLM WITH SAVING SPECIES INFLUX FILE");
	}
        if(tmpEnvironment->getEnergy() <= 1)
        {
            if(!tmpEnvironment->saveNrgBoolFncStructureSTD(tmpSavingPath))
                    ExitWithError("saveNrgBoolFncStructure", "PROBLM WITH SAVING ENERGY BOOLEAN FUNCTIONS FILE");
        }
    if(!tmpEnvironment->saveSpeciesStructureSTD(tmpGen, tmpSim, tmpStep, tmpSavingPath))
		ExitWithError("saveSpeciesStructure", "PROBLEM WITH SAVING INITIAL SPECIES STRUCTURE TO FILE");
    if(!tmpEnvironment->saveReactionsStructureSTD(tmpGen, tmpSim, tmpStep, tmpSavingPath))
		ExitWithError("saveReactionsStructure", "PROBLEM WITH SAVING INITIAL REACTIONS STRUCTURE TO FILE");
    if(!tmpEnvironment->saveCatalysisStructureSTD(tmpGen, tmpSim, tmpStep, tmpSavingPath))
		ExitWithError("saveCatalysisStructure", "PROBLEM WITH SAVING INITIAL CATALYSIS STRUCTURE TO FILE");
}
