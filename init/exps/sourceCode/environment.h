/** \class environment
 * \brief environment class
 *
 * This class contains environmental proprieties and methods; within this class all things occur
 * \author Alessandro Filisetti
 * \version 0.2
 * \date 2011-12-15
 */

#ifndef ENVIRONMENT
#define ENVIRONMENT

#include "acs_headers.h"
#include "species.h"
#include "reactions.h"
#include "commonFunctions.h"
#include "catalysis.h"
#include "gillespie.h"

class environment
{
private:
    acs_int nSIM; // Number of simulations
    acs_int nGEN; // Number of generations
    acs_double actualTime; // Total Simulation Real Time
    acs_double nSeconds; // Number of physical seconds
    acs_int nReactions; // MAX Number of reactions
    acs_double nHours; // MAX number of hours
    acs_int nAttempts; // MAX number of attempts for simulation
    acs_int Currentattempts; // Current number of temptatives
    acs_double timeStructuresSavingInterval; // Time between two successive files storing procedure
    acs_double fileTimesSaveInterval; // Time between two successive times.csv storing procedure
    acs_double newSpeciesProbMinThreshold; // Minimal new species creation probability to allow system expansion
    acs_int lastFiringDiskSpeciesID; // ID of the last species of the firing disk
    acs_double overallConcentration; // initial overall concentration
    acs_double ECConcentration; // Energie Carriers concentration
    acs_int nonCatalyticMaxLength; // Max length of the NON catalytic species
    acs_longInt nonCatalyticLastID; // ID of the non catalytic last species.
    acs_double reactionProbability; // Reaction probability
    acs_double cleavageProbability; // Cleavage probability (when new reactions are created)
    bool main_rev_rct_allowed; // Possibility to have both forward and backward reactions as main reaction
    bool reverseReactions; // if true backward reactions are present
    acs_double revRctRatio; // ratio between forward and backward reactions, it is used only in the creation of new reactions (if reverseReactions = TRUE)
    acs_double spontRct; // if true spontaneous reactions are present
    species nrgCarrier; // !!! declare energy carrier, !!! PRESENTLY this is not used... it was used in previous versions !!!
    acs_longInt numberOfMolecules; // Total number of molecules
    acs_longInt numberOfNewMolecules; // Total number of NEW molecules
    acs_longInt numberOfSpecies; // Total number of species
    acs_longInt numberOfNewSpecies; // Total number of NEW species
    acs_longInt numberOfCpx; // Total number of complexes
    acs_longInt numberOfCpxMols; // Total number of new complexes
    acs_longInt cleavageCounter;            // Number of Cleavage counter
    acs_longInt endoCleavageCounter;        // Number of endoergonic cleavage counter
    acs_longInt condensationCounter;        // Number of condensation counter
    acs_longInt endoCondensationCounter;   // Number of endoergonic condensation counter
    acs_longInt cpxFormCounter;	   // Complex formation Counter
    acs_longInt cpxDissCounter;		// Complex dissociation counter
    acs_longInt spontDissCounter; // Spontaneous dissociations counter
    acs_longInt spontAssCounter; // Spontaneous Associations counter
    int overallLoadedMolsCounter;  // Total number of loaded molecules
    acs_double internalTimesStoredCounter; // Time interval to save certain structures
    acs_int lastEvaluatedSpeceisForNewReactions; // Last species evaluated for new reactions, in such a way once that a species has been evaluated to be involved in new reactions

    // Kinetic Constants --------------------------

    acs_double K_ass;
    acs_double K_diss;
    acs_double K_cpx;
    acs_double K_cpxDiss;
    acs_double K_nrg;
    acs_double K_nrg_decay; // Energy Decay (both for species and carriers) kinetic constant
    acs_double K_spont_diss; // Spontaneous  dissociation kinetic constant
    acs_double K_spont_ass; // Spontaneous association kinetic constant
    acs_int    maxLOut; // Maximum length of the species that may be involved in the outflow process
    acs_double cleavage_KineticConstant;
    acs_double complex_KineticConstant;
    acs_double condensation_KineticConstant;
    acs_double complexDeg_KineticConstant;
    acs_double moleculeDecay_KineticConstant;

    // SYSTEM PARAMETERS --------------------------
    acs_int energy; // Energy control variable
    acs_double ratioSpeciesEnergizable; // Proabability for a species to be energizable.
    acs_int nrgBoolFlag; // variable indicating the presence (1) or not (0) of the energy
    acs_int solubility_threshold;
    acs_double diffusion_contribute;
    acs_double influx_rate; // how many molecules per second are introduced in the system
    acs_double timeSinceTheLastInFlux; // Time elapsed from the previous injection of nutrients
    string alphabet; // Alphabet (e.g. AB, RNA, DNA, Proteins)
    acs_double volume; // Reactor volume
    bool volumeGrowth; // Boolean variable defining the possibility for the volume to change
    bool stochDivision; // Boolean Variable defining the division process (stochastic or deterministic)
    acs_double randomSeed; //Random seed of the simulation
    //TR acs_double totalPossibleNumberOfSpecies; // Overall number of species in the firing disk. According to this parameters the initial set of molecules is set
    vector<species> allSpecies; // List containing all the species
    vector<species> firingDisk; // Species present in the firing disk
    vector<acs_int> nutrientsForInflux; // Create an influx species IDs list
    vector<acs_double> nutrientsProb2BeSelected; // create an influx parallel list containing the probability to be selected for the influx species
    vector<reactions> allReactions; // List containing all the reactions
    vector<catalysis> allCatalysis; // List containing all the catalysis
    vector<species> initialAllSpecies; // initial species population when generation is graeater than 1
    vector<reactions> initiaAllReactions; // initial reactions when generation is graeater than 1
    vector<catalysis> initialAllCatalysis; // initial catalysis when generation is graeater than 1
    vector<gillespie> allGillespieScores; // List containing the Gilespie structure
    vector<gillespie> COPYOFallGillespieScores; // List containing a copy of Gilespie structure to save to file
    vector<acs_double> gillespieCumulativeStepScoreList; // List containing cumulative gillespie score values
    vector<acs_double> reactionsTime; // List of all reaction time
    vector<acs_longInt> gillespieReactionsSelected; // List of the reactions selected by means of the Gillespie algorithm
    vector<acs_double> allTimes; // List of the computational times of each reaction
    int debugLevel; // level of message details show during the simution
    vector<acs_int> livingSpeciesIDsList; // this list contains all the living species ID. It will be updated after each efflux action
    vector<acs_double> cumLivingAmountsList; // this list contains the cumulative number of molecules according to the LivingSpeciesIDslist
    vector<acs_int> unchargedIDlist; // This list contains all the uncharged species ID
    vector<acs_double> cumUnchargedAmountList; // This list contains the cumulative numbers of uncharged molecule. Thi list is correlated with the unchargedIdList
    vector<acs_int> chargedIDlist; // This list contains all the charged speciesID
    vector<acs_double> cumChargedAmountList; // This List contains all the charged speciesID
    vector<acs_int> nrgBooleanFunctions; // This list contains all the possible energetic Boolean Functions
    vector<acs_double> nrgBoolFncsProb2BeSelected; // create a Boolean functions parallel list containing the probability for a Bool fnc to be selected
	
    // SUPPORT VARIABLES ----------------------------
	//TR acs_double decimalCarrierToLoad; // Decimal part of carrier to load in the next reaction
	acs_double decimalMoleculesToEfflux; // Decimal part of molecule to efflux in the next reaction
	acs_double decimalMoleculesToLoad; // Decimal part of molecules to load in the next reaction
	acs_double decimalMoleculesToUNLOAD; // Decimal part of molecules to unload because of the energy decay 
	acs_double decimalComplexesToDissociate; // Decimal part of complex to dissociate.
	
    // STATISTICAL VARIABLES ------------------------
	acs_double gillespieMean;       // Mean of all the Gillespie Scores
    acs_double gillespieSD;         // Standard deviation of all the Gillespie Scores
    acs_double gillespieEntropy;    // sum of all the p(log(p)) of the Gillespie Scores
	acs_double gillespieTotalScore; // Gillespie total score to reset at each gillespie computation
    acs_double gillespieNewSpeciesScore; // The sum of the gillespie scores creating only new species
    acs_double ratioBetweenNewGillTotGill; // gillespieNewSpeciesScore / gillespieTotalScore
    acs_double reverseReactionsGillScore; // the sum of the gillespie scores related to the reverse reactions
    acs_double ratioBetweenReverseAndTotalScore; // Ratio between the reverse reaction scores and all the reaction scores
    clock_t gillespiePartialTimer; // Gillespie algorithm computational time
    vector<acs_double> gillespiePartialTimes;
    clock_t performReactionPartialTimer; // The time necessary for the evaluation of a new species
    vector<acs_double> performReactionPartialTimes;
    clock_t remainingProcessesPartialTimer; // THe time needed for all the other procedure
    vector<acs_double> remainingProcessesPartialTimes;

	// METHODS ---------------------------|
	
	public:
	//environment();
    environment(string tmpInitialPath);
	~environment(){}
	
	//GET FUNCTIONS
	acs_int getNgen()const{return nGEN;}
	acs_int getNsim()const{return nSIM;}
	acs_double getActualTime()const{return actualTime;}
	acs_double getNseconds()const{return nSeconds;}
	acs_int getNreactions()const{return nReactions;}
    acs_double getMAXhours()const{return nHours;}
    acs_int getMAXattempts()const{return nAttempts;}
    acs_int getCurrentAttempts()const{return Currentattempts;}
	acs_double getTimeStructuresSavingInterval()const{return timeStructuresSavingInterval;}
    acs_double getFileTimesSavingInterval()const{return fileTimesSaveInterval;}
	acs_int getLastFiringDiskSpeciesID()const{return lastFiringDiskSpeciesID;}
	acs_double getOverallConcentration()const{return overallConcentration;}
	acs_int getMaxNonCatalyticLength()const{return nonCatalyticMaxLength;}
	acs_double getRctProb()const{return reactionProbability;}
	acs_double getCleavProb()const{return cleavageProbability;}
	bool getReverseReactions()const{return reverseReactions;}
	acs_int getEnergy()const{return energy;}
	acs_double getRatioSpeciesEnergizable()const{return ratioSpeciesEnergizable;}
	acs_int getADP()const{return nrgCarrier.getNOTchargeMols();}
	acs_int getATP()const{return nrgCarrier.getChargeMols();}
	acs_longInt getMols()const{return numberOfMolecules;}
	acs_longInt getNewMols()const{return numberOfNewMolecules;}
	acs_longInt getNspecies()const{return numberOfSpecies;}
	acs_longInt getNnewSpecies()const{return numberOfNewSpecies;}
	acs_longInt getNcpx()const{return numberOfCpx;}
	acs_longInt getNcpxMols()const{return numberOfCpxMols;}
	
	acs_double getGillespieMean()const{return gillespieMean;}
	acs_double getgillespieSD()const{return gillespieSD;}
	acs_double getgillespieEntropy()const{return gillespieEntropy;}
    acs_double getRatioBetweenNewGillTotGill()const{return ratioBetweenNewGillTotGill;}
    acs_double getRatioBetweenBackandForw()const{return ratioBetweenReverseAndTotalScore;}
    bool getSystemExpFlag()const{return newSpeciesProbMinThreshold < ratioBetweenReverseAndTotalScore;}
	
	// Kinetic Constants -----------------------------
	acs_double getKdiss()const{return K_diss;}
	acs_double getKass()const{return K_ass;}
	acs_double getKcpx()const{return K_cpx;}
	acs_double getKcpxDiss()const{return K_cpxDiss;}
	acs_double getKnrg()const{return K_nrg;}
	acs_double getKirrad()const{return K_nrg_decay;}
	acs_double getK_spont_diss()const{return K_spont_diss;}
	acs_double getK_spont_ass()const{return K_spont_ass;}
	
	acs_double getCleavageKC()const{return cleavage_KineticConstant;}
	acs_double getComplexKC()const{return complex_KineticConstant;}
	acs_double getCondensationKC()const{return condensation_KineticConstant;}
	acs_double getComplexDegKC()const{return complexDeg_KineticConstant;}
	acs_double getMoleculeDecayKC()const{return moleculeDecay_KineticConstant;}
	acs_int getMaxLOut()const{return maxLOut;}
	acs_int getSolubilityThreshold()const{return solubility_threshold;}

	// General Getting function
	acs_double getDiffusionContribute()const{return diffusion_contribute;}
	acs_double getInflux()const{return influx_rate;}
	acs_double getRefillInterval()const{return influx_rate;}
	string getAlphabet()const{return alphabet;}
	acs_double getVolume()const{return volume;}
	acs_double getRandomSeed()const{return randomSeed;}
    vector<species> getMoleculesPopulation()const{return allSpecies;}
	acs_longInt getTotalNumberOfSpecies();
	acs_longInt getTotalNumberOfPossibleCatalysts();
	acs_longInt getTotalNumberOfMolecules();
	acs_longInt getTotalNumberOfComplexSpecies();
	acs_longInt getTotalNumberOfComplexes();
	acs_longInt getTotalNumberOfMonomers();
    vector<reactions> getReactionsLayer()const{return allReactions;}
	int getDebugLevel()const{return debugLevel;}
	acs_longInt getNumberOfTheoreticalSpecies()const{return (acs_longInt)allSpecies.size();}
	acs_longInt getNumberOfReactions()const{return (acs_longInt)allReactions.size();}
	acs_longInt getNumberOfCatalysis()const{return (acs_longInt)allCatalysis.size();}
	acs_longInt getNumberOfGillespieCOPYpossibleRcts()const{return (acs_longInt)COPYOFallGillespieScores.size();}
	acs_longInt getNumberOfGillespiePossibleRcts()const{return (acs_longInt)allGillespieScores.size();}
	void setLivingSpeciesIDsAndAmounts();
	void setNotChargedAndChargedSpeciesIDsAndAmounts();
	acs_longInt getCleavageCounter()const{return cleavageCounter;}
	acs_longInt getEndoCleavageCounter()const{return endoCleavageCounter;}
	acs_longInt getCondensationCounter()const{return condensationCounter;}
	acs_longInt getEndoCondensationCounter()const{return endoCondensationCounter;}
	acs_longInt getCpxFormCounter()const{return cpxFormCounter;}
	acs_longInt getCpxDissCounter()const{return cpxDissCounter;}
	acs_longInt getOverallLoadedMolsCounter()const{return overallLoadedMolsCounter;}
	acs_longInt getSpontDissCounter()const{return spontDissCounter;}
	acs_longInt getSpontAssCounter()const{return spontAssCounter;}
	acs_int getTotNumberOfChargedMols();

	// PROMPT FUNCTIONS
	void showGlobalParameter();
	void printInitialCondition();
	void printAllSpeciesIdAndSequence();
	void printGillespieStructure();
	void printNutrientsAndProbability();
	void printAllChargeMols();
	
    // FILE STRUCTURE UPLOAD FUNCTION
    bool createInitialMoleculesPopulationFromFileSTD(string tmpSpeciesFilePath);
    bool createInitialReactionsLayerFromFileSTD(string tmpSpeciesFilePath);
    bool createInitialCatalysisLayerFromFileSTD(string tmpCatalysisFilePath);
    bool createInfluxLayersFromFileSTD(string tmpInfluxFilePath);
    bool createNrgBooleanFunctionsFromFileSTD(string tmpBoolNrgFilePath);

    bool createInitialMoleculesPopulationFromSpecificFileSTD(string tmpSpeciesFilePath, acs_int tmpActGEN, acs_int tmpActSIM);
    bool createInitialReactionsLayerFromSpecificFileSTD(string tmpReactionsFilePath, acs_int tmpActGEN, acs_int tmpActSIM);
    bool createInitialCatalysisLayerFromSpecificFileSTD(string tmpCatalysisFilePath, acs_int tmpActGEN, acs_int tmpActSIM);
	
	// OPERATIVE FUNCTIONS
	//TR void createInitialMoleculesPopulation(MTRand& tmpRndDoubleGen);
	//TR void createReactionsLayer(MTRand& tmpRndDoubleGen);
	void nutrientsAmountsFixing();
	acs_int computeSngSpeciesRctsNumber(acs_longInt tmpTotalNumberOfReactions, MTRand& tmpRndDoubleGen);
	acs_int selectWhetherCleavageOrCond(MTRand& tmp__RndDoubleGen);
    bool createReactionsForThisSpecies(acs_longInt tmpsID, acs_int tmpReactionsForThisSpecies, MTRand& tmp_RndDoubleGen, vector<acs_longInt>& tmpIDOfCandidateSpecies, acs_int tmpRctCreationType);
    bool updateReactions(acs_longInt tmpIDtoUpdate, acs_longInt tmpNewSpecies, acs_int tmpRctType, vector<acs_longInt>& tmp_AlreadyEvaluatedSpeciesVector,  MTRand& tmp_RndDoubleGen);
	acs_double createDiffusionRenforcement(acs_double tmpDiffEnh, acs_int tmpNewSpeciesLength);
	bool setSolubility(acs_int tmpNewSpeciesLength, MTRand& tmpRndDoubleGen);
	//TR acs_longInt createInitialAmount(int tmpAlphabetLength, int tmpSpeciesLength);
	//TR acs_double createInitialConcentration(int tmpAlphabetLength, int tmpSpeciesLength);
    acs_longInt returnPosSpeciesAlreadyPresent(string tmpNewSequence); //return vector species size if it doesn't
    acs_longInt returnPosReactionAlreadyPresent(acs_int tmpReactionType, acs_longInt tmpIds_I, acs_longInt tmpIds_II, acs_longInt tmpIds_III); //return vector allReactions size if it doesn't
	bool checkIfTheReactionIsAlreadyCatalyzedByThisSpecies(acs_longInt tmpSPeciesID, acs_longInt tmpIdReaction);
    bool performOPTGillespieComputation(MTRand& tmpRndDoubleGen, clock_t& tmpTimeElapsed, acs_int tmpActGEN, acs_int tmpActSIM, acs_int tmpActSTEP, string tmpStoringPath);
    bool performReaction(acs_longInt reaction_u, MTRand& tmp_RndDoubleGen, acs_int tmp_ActGEN, acs_int tmp_ActSIM, acs_int tmp_ActSTEP, string tmp_StoringPath);
	bool newSpeciesEvaluationIII(acs_longInt tmpNewSpecies, MTRand& tmp___RndDoubleGen);
    bool complexEvaluation(string tmpComplex, MTRand& tmp___RndDoubleGen, acs_int tmpCuttingPnt, acs_longInt tmpCatalyst_ID, acs_longInt tmpCatID, acs_longInt tmpSubstrate_ID, acs_longInt tmpSecSub_ID, bool tmpCpxType);
	acs_double computeSinglGilScore(acs_longInt tmpAmountI, acs_double tmpDifI, acs_int tmpSolI,acs_longInt tmpAmountII, acs_double tmpDifII, acs_int tmpSolII, acs_double tmpK, bool tmpSameMol);
	void performSingleGilleSpieIntroduction(acs_longInt tmpAmountI, acs_longInt tmpAmountII, acs_longInt tmpIDI, acs_longInt tmpIDII, acs_longInt tmpIDCatalysis, acs_int tmp__rctType,
											acs_longInt tmpMol_I, acs_longInt tmpMol_II, acs_longInt tmpMol_III, acs_longInt tmpMol_IV, acs_int tmpNRGDirection, acs_longInt tmpRctID,
											bool tmpSameSpeciesControl);
	// VOLUME OPERATIVE FUNCTIONS
	void changeVolume(acs_int tmpTimeSinceLastReaction);
	
	// Change overall COUNTER number of species, molecules, complexes type and complexes token
	void incNumberOfSpecies(){numberOfSpecies++;}
	void decNumberOfSpecies(acs_int tmpID){if(allSpecies.at(tmpID).getAmount() == 0) numberOfSpecies--; }
	void incNumberOfMols(){numberOfMolecules++;}
	void decNumberOfMols(){numberOfMolecules--;}
	void incNumberOfCpx(){numberOfCpx++;}
	void decNumberOfCpx(acs_int tmpID){if(allSpecies.at(tmpID).getAmount() == 0) numberOfCpx--;}
	void incNumberOfCpxMols(){numberOfCpxMols++;}
	void decNumberOfCpxMols(){numberOfCpxMols--;}
	
	void incNumberOfNewSpecies(acs_int tmpID){if(tmpID > lastFiringDiskSpeciesID) numberOfNewSpecies++;}
	void decNumberOfNewSpecies(acs_int tmpID){if((tmpID > lastFiringDiskSpeciesID)&&(allSpecies.at(tmpID).getAmount() == 0)) numberOfNewSpecies--;}
	void incNumberOfNewMols(acs_int tmpID){if(tmpID > lastFiringDiskSpeciesID) numberOfNewMolecules++;}
	void decNumberOfNewMols(acs_int tmpID){if(tmpID > lastFiringDiskSpeciesID) numberOfNewMolecules--;}
	
	void decMolSpeciesProcedure(acs_int tmp_ID){decNumberOfSpecies(tmp_ID); decNumberOfMols(); decNumberOfNewSpecies(tmp_ID); decNumberOfNewMols(tmp_ID);}
	void decCpxProcedure(acs_int tmp_ID){decNumberOfCpx(tmp_ID);decNumberOfCpxMols();}
	void incMolProcedure(acs_int tmp_ID){incNumberOfMols();incNumberOfNewMols(tmp_ID);}
	void incSpeciesProcedure(acs_int tmp_ID){if(allSpecies.at(tmp_ID).getAmount() == 1){incNumberOfNewSpecies(tmp_ID); incNumberOfSpecies();}}
	//void chargeMolProcess(acs_int tmp_ID){}
    void unchargeMolProcess(acs_int tmp_ID){allSpecies[tmp_ID].unchargeMol();decOverallLoadedMolsCounter();}// removeChargeMolFromList(tmp_ID);}
	
	void incCleavageCounter(){cleavageCounter++;}
	void incEndoCleavageCounter(){endoCleavageCounter++;}
	void incCondensationCounter(){condensationCounter++;}
	void incEndoCondensationCounter(){endoCondensationCounter++;}
	void incCpxFormCounter(){cpxFormCounter++;}
	void incCpxDissCounter(){cpxDissCounter++;}
	void incOverallLoadedMolsCounter(){overallLoadedMolsCounter++;}
	void decOverallLoadedMolsCounter(){overallLoadedMolsCounter--;}
	void incSpontDissCounter(){spontDissCounter++;}
	void incSpontAssCounter(){spontAssCounter++;}

	void resetCleavageCounter(){cleavageCounter = 0;}
	void resetEndoCleavageCounter(){endoCleavageCounter = 0;}
	void resetCondensationCounter(){condensationCounter = 0;}
	void resetEndoCondensationCounter(){endoCondensationCounter = 0;}
	void resetOverallLoadedMolsCounter(){overallLoadedMolsCounter = 0;}
	void resetCpxFormCounter(){cpxFormCounter = 0;}
	void resetCpxDissCounter(){cpxDissCounter = 0;}
	void resetSpontDissCounter(){spontDissCounter = 0;}
	void resetSpontAssCounter(){spontAssCounter = 0;}

	void resetReactionsCounter(){resetCleavageCounter(); resetEndoCleavageCounter();
		resetCondensationCounter(); resetEndoCondensationCounter(); resetOverallLoadedMolsCounter();
		resetCpxFormCounter(); resetCpxDissCounter(); resetSpontDissCounter(); resetSpontAssCounter();}

	bool addChargeMolToList(acs_int tmpSpeciesID);
	bool removeChargeMolFromList(acs_int tmpSpeciesID);

    // STRUCTURE CLEANING
	void clearAllStructures();
    //void resetConcentrationToInitialConditions(QString tmpSpeciesFilePath, acs_int tmpActGEN, acs_int tmpActSIM);
    void resetConcentrationToInitialConditions();
    void storeInitialStructures();

    //  DETERMINISTIC FUNCTIONS
	bool performRefill(acs_double tmpTimeSinceTheLastInFlux, acs_double tmpMinimalTimeForOneMols, MTRand& tmp__RndDoubleGen);
	bool performMoleculesEfflux(acs_double tmpTimeInterval, MTRand& tmp_RndDoubleGen);
	bool performDETMoleculesCharging(acs_double tmpTimeInterval, MTRand& tmp_RndDoubleGen);
	bool performDETComplexDissociation(acs_double tmpTimeInterval, MTRand& tmp_RndDoubleGen);
	void setActualTime(acs_double tmpActualTime){actualTime = tmpActualTime;}
	//TR void getNutrientsFromTheFiringDisk();
	void updateSpeciesAges();
    void increaseAttempts(){Currentattempts++;}
	
	// REACTION PERFORMING FUNCTIONS 
	bool performCondensation(acs_longInt tmpCatalyst, acs_longInt tmpSubstrate, acs_longInt tmpProduct, acs_longInt tmpComplex,
							 acs_longInt tmpIdReaction, acs_longInt tmpIdCatalysis, MTRand& tmp__RndDoubleGen);
	bool perform_endo_Condensation(acs_longInt tmpCatalyst, acs_longInt tmpSubstrate, acs_longInt tmpProduct, acs_longInt tmpComplex,
								   acs_int tmpNRGside, acs_longInt tmpIdReaction, acs_longInt tmpIdCatalysis, MTRand& tmp__RndDoubleGen);
	bool performCleavage(acs_longInt tmpSubstrate, acs_longInt tmpProduct_I, acs_longInt tmpProduct_II, 
						 acs_longInt tmpIdReaction, acs_longInt tmpIdCatalysis, MTRand& tmp__RndDoubleGen);
	bool perform_endo_Cleavage(acs_longInt tmpSubstrate, acs_longInt tmpProduct_I, acs_longInt tmpProduct_II,
							   acs_int tmpNRGside, acs_longInt tmpIdReaction, acs_longInt tmpIdCatalysis, MTRand& tmp__RndDoubleGen);
    bool performComplexFormation(acs_longInt tmpCatalyst, acs_longInt tmpSubstrate, acs_longInt tmpCatID, acs_longInt tmpSecSub, MTRand& tmp__RndDoubleGen);
    bool perform_endo_ComplexFormation(acs_longInt tmpCatalyst, acs_longInt tmpSubstrate, acs_longInt tmpCatID, acs_longInt tmpSecSub, acs_int tmpNRGSide, MTRand& tmp__RndDoubleGen);
	bool performComplexDissociation(acs_longInt tmpComplex, acs_longInt tmpCatalyst, acs_longInt tmpSubstrate, MTRand& tmp__RndDoubleGen);
	bool performSpontaneousCondensation(acs_longInt tmpReaction, MTRand& tmp__RndDoubleGen);
	bool performSpontaneousCleavage(acs_longInt tmpReaction, MTRand& tmp__RndDoubleGen);
	bool performMoleculeEfflux(acs_longInt tmpSpecies, MTRand& tmp__RndDoubleGen);
    //TR bool performPhosphorilation(acs_longInt tmpSpecies);
	bool performEnergyEfflux(MTRand& tmp__RndDoubleGen);
	
	// COHERENCE CONTROL FUNCTION
	bool structureCoherenceCheckUp();
	bool notInverseReactionAlreadyCatalyzed(acs_int tmpRct, acs_longInt tmpID_I, acs_longInt tmpID_II);
	bool checkIfOnlyMutualCatalysis(acs_int tmpCat, acs_int tmpCandidateProduct);
    bool checkAvailability(acs_longInt tmpMI, acs_longInt tmpMII, acs_longInt tmpQI, acs_longInt tmpQII);
    void inserSubListInSpecies();
    void showSubListInSpecies();
    void showGillEngagementInSpecies();
	
    // SAVE TO FILE FUNCTIONS

    bool saveConfigurationFileSTD(string tmpStoringPath);
    bool saveInfluxStructureSTD(string tmpStoringPath);
    bool saveNrgBoolFncStructureSTD(string tmpStoringPath);
    string zeroBeforeStringNumberSTD(acs_int tmpTotN, acs_int tmpCurrentN);
    bool saveSpeciesStructureSTD(acs_int tmpCurrentGen, acs_int tmpCurrentSim, acs_int tmpCurrentStep, string tmpStoringPath);
    bool saveReactionsStructureSTD(acs_int tmpCurrentGen, acs_int tmpCurrentSim, acs_int tmpCurrentStep, string tmpStoringPath);
    bool saveCatalysisStructureSTD(acs_int tmpCurrentGen, acs_int tmpCurrentSim, acs_int tmpCurrentStep, string tmpStoringPath);
    bool saveTimesSTD(acs_int tmpCurrentGen, acs_int tmpCurrentSim, acs_int tmpCurrentStep, string tmpStoringPath);
    bool saveReactionsParametersSTD(acs_int tmp__CurrentGen, acs_int tmp__CurrentSim, acs_int tmp__CurrentStep, string tmp__StoringPath, acs_int tmpRctType,
                                 acs_longInt tmpCat, acs_longInt tmpMol_I, acs_longInt tmpMol_II, acs_longInt tmpMol_III);
    bool saveLivingSpeciesIDSTD(acs_int tmp__CurrentGen, acs_int tmp__CurrentSim, acs_int tmp__CurrentStep, string tmp__StoringPath);
    bool saveLivingSpeciesAmountSTD(acs_int tmp__CurrentGen, acs_int tmp__CurrentSim, string tmp__StoringPath);
    bool saveLivingSpeciesConcentrationSTD(acs_int tmp__CurrentGen, acs_int tmp__CurrentSim, string tmp__StoringPath);
	//bool saveGillespieStructure(acs_int tmpCurrentSim, acs_int tmpCurrentStep, QString tmpStoringPath);
	//bool saveAdjacentStructures(acs_int tmpCurrentSim, acs_int tmpCurrentStep, QString tmpStoringPath);

    // STATISTIC FUNCTIONS
    bool devStd();
    bool entropy();

};


#endif
