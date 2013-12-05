/**  \class environment
 * \author Alessandro Filisetti
 * \version 2.4
 * \date 2010-06-10
 * *^*^*^*^*
 */
#include "environment.h"

/**
 * TEMPLATE DEACLARATION AND DEFINITION
 */

/* ************
|
|	CONSTRUCTORS
|
****************/

/**
 Test environment costructor
 @version 2.4
 @date 2010-06-27
 @param tmpRndDoubleGen randomGenerator reference
*/
/*environment::environment()
{
	//TR onlyEnvironmentCreation = 0;
	nGEN = 1;
	nSIM = 1;
	nReactions = 100;
	//TR initialPopulationNumber = 30;
	timeStructuresSavingInterval = 1;
   	fileTimesSaveInterval = 0;
   	newSpeciesProbMinThreshold = 0;
	//TR initialMaxLength = 4;
	reactionProbability = 0.01;
	cleavageProbability = 0.5;
	cleavage_KineticConstant = 1.5;
	complex_KineticConstant = 0.01;
	condensation_KineticConstant = 0.5;
	complexDeg_KineticConstant = 0.001;
	moleculeDecay_KineticConstant = 0.1;
	solubility_threshold = 0;
	diffusion_contribute = 1;
	//TR influx = 1;
	influx_rate = 0.5;
	alphabet = "AB";
	randomSeed = 100;
	debugLevel = 1;
	timeSinceTheLastInFlux = 0;
	actualTime = 0;
	energy = 0;
	ratioSpeciesEnergizable = 0;
	// --------------------------
	decimalMoleculesToLoad = 0;
	decimalMoleculesToUNLOAD = 0;
	decimalMoleculesToEfflux = 0;
	decimalComplexesToDissociate = 0;
	numberOfSpecies = 0;
	numberOfNewSpecies = 0;
	numberOfMolecules = 0;
	numberOfNewMolecules = 0;
	numberOfCpx = 0;
	numberOfCpxMols = 0;
	gillespieTotalScore = 0;
	maxLOut = 0;
    internalTimesStoredCounter = 0;
    nHours = 5;
    nAttempts = 5;
    Currentattempts = 0;
    revRctRatio = 1000;
}*/

/**
 Environment Constructor
 @version 1.0
 @param tmpInitialPath
 @date 2013/07/04
 */
environment::environment(string tmpInitialPath)
{
    if(debugLevel == FINDERRORDURINGRUNTIME) cout << "environment::environment start" << endl;

    // CONFIGURATION FILE PATH CREATION
    string confFilePath = tmpInitialPath + "acsm2s.conf";
    main_rev_rct_allowed = 0; //Default value if not present in old conf files
    ifstream myfile;
    myfile.open(confFilePath.c_str());
    string strID;
    vector<string> linered;
    while (myfile.good())
    {
        getline(myfile, strID);
        if(!strID.empty())
            if(strID.c_str()[0] != '#'){
                //cout << strID << endl;
                linered = split(strID, "=");
                if(linered[0] == "nGEN") nGEN = atoi(linered[1].c_str());
                if(linered[0] == "nSIM") nSIM = atoi(linered[1].c_str());
                if(linered[0] == "nSeconds") nSeconds = atof(linered[1].c_str());
                if(linered[0] == "nReactions") nReactions = atoi(linered[1].c_str());
                if(linered[0] == "randomSeed") randomSeed = atoi(linered[1].c_str());
                if(linered[0] == "debugLevel") debugLevel = atoi(linered[1].c_str());
                if(linered[0] == "timeStructuresSavingInterval") timeStructuresSavingInterval = atof(linered[1].c_str());
                if(linered[0] == "fileTimesSaveInterval") fileTimesSaveInterval = atof(linered[1].c_str());
                if(linered[0] == "newSpeciesProbMinThreshold") newSpeciesProbMinThreshold = atof(linered[1].c_str());
                if(linered[0] == "nHours") nHours = atof(linered[1].c_str());
                if(linered[0] == "nAttempts") nAttempts = atoi(linered[1].c_str());

                // ENVIRONMENTAL PARAMETERS
                if(linered[0] == "lastFiringDiskSpeciesID") lastFiringDiskSpeciesID = atoi(linered[1].c_str()); // DA ELIMINARE CON ATTENZIONE
                if(linered[0] == "overallConcentration") overallConcentration = atof(linered[1].c_str()); // DA ELIMINARE CON ATTENZIONE
                if(linered[0] == "ECConcentration") ECConcentration = atof(linered[1].c_str()); // DA TRASFORMARE NELLA VARIABILE INDICANTE LA QUANTITA DI CARRIER SEMPRE PRESENTI
                if(linered[0] == "alphabet") alphabet = linered[1];
                if(linered[0] == "volume") volume = atof(linered[1].c_str());
                if(linered[0] == "volumeGrowth") volumeGrowth = fromStrToBool(linered[1]);
                if(linered[0] == "stochDivision") stochDivision = fromStrToBool(linered[1]);

                // DYNAMIC PARAMETERS
                if(linered[0] == "energy") energy = atoi(linered[1].c_str()); // DA ELIMINARE CON ATTENZIONE, IL TUTTO VERRA' SOSTITUITO DALLA FUNZIONE BOOLEANA
                if(linered[0] == "ratioSpeciesEnergizable") ratioSpeciesEnergizable = atof(linered[1].c_str());
                if(linered[0] == "nonCatalyticMaxLength") nonCatalyticMaxLength = atoi(linered[1].c_str());
                if(linered[0] == "reactionProbability") reactionProbability = atof(linered[1].c_str());
                if(linered[0] == "cleavageProbability") cleavageProbability = atof(linered[1].c_str());
                if(linered[0] == "main_rev_rct_allowed") main_rev_rct_allowed = atoi(linered[1].c_str());
                if(linered[0] == "reverseReactions") reverseReactions = atoi(linered[1].c_str());
                if(linered[0] == "revRctRatio") revRctRatio = atof(linered[1].c_str());
                if(linered[0] == "spontRct") spontRct = atof(linered[1].c_str());
                if(linered[0] == "K_ass") K_ass = atof(linered[1].c_str());
                if(linered[0] == "K_diss") K_diss = atof(linered[1].c_str());
                if(linered[0] == "K_cpx") K_cpx = atof(linered[1].c_str());
                if(linered[0] == "K_cpxDiss") K_cpxDiss = atof(linered[1].c_str());
                if(linered[0] == "K_nrg") K_nrg = atof(linered[1].c_str()); // DA ELIMINARE CON ATTENZIONE
                if(linered[0] == "K_nrg_decay") K_nrg_decay = atof(linered[1].c_str()); // DA ELIMINARE CON ATTENZIONE
                if(linered[0] == "moleculeDecay_KineticConstant") moleculeDecay_KineticConstant = atof(linered[1].c_str());
                if(linered[0] == "K_spont_ass") K_spont_ass = atof(linered[1].c_str());
                if(linered[0] == "K_spont_diss") K_spont_diss = atof(linered[1].c_str());
                if(linered[0] == "maxLOut") maxLOut = atoi(linered[1].c_str());
                if(linered[0] == "solubility_threshold") solubility_threshold = atoi(linered[1].c_str());
                if(linered[0] == "diffusion_contribute") diffusion_contribute = atof(linered[1].c_str());
                if(linered[0] == "influx_rate") influx_rate = atof(linered[1].c_str());
            }
    }
    myfile.close();
    //OPEN FILE

    if(nSeconds < timeStructuresSavingInterval)
        ExitWithError("environment::environment","No Structures will be saved during the simulation");

    if(nSeconds < fileTimesSaveInterval)
        ExitWithError("environment::environment","No times file will be saved during the simulation");

    if(debugLevel >= RUNNING_VERSION)
        showGlobalParameter();
    //IF RANDOM SEED IS 0 IT HAS TO BE CREATED RANDOMLY
    if(randomSeed == 0)
    {
        //randomSeed = getIntRandom(1, RANDOMRANGE, tmpRndDoubleGen);
        srand(time(0));
        randomSeed = rand();
        if(debugLevel >= RUNNING_VERSION)
            cout << "|- New random seed -> " << randomSeed << endl;
    }

    if(energy == 1)
        nrgBoolFlag = ENERGYBASED;
    else
        nrgBoolFlag = ENERGYFREE;

    // INITIALIZE THE REFILL TIMER
    timeSinceTheLastInFlux = 0;
    // INITIALIZE TIME TO 0
    actualTime = 0;

    numberOfSpecies = 0;
    numberOfNewSpecies = 0;
    numberOfMolecules = 0;
    numberOfNewMolecules = 0;
    numberOfCpx = 0;
    numberOfCpxMols = 0;
    decimalMoleculesToEfflux = 0;
    decimalMoleculesToLoad = 0;
    decimalMoleculesToUNLOAD = 0;
    decimalComplexesToDissociate = 0;
    gillespieTotalScore = 0;
    gillespieNewSpeciesScore = 0;
    ratioBetweenNewGillTotGill = 0;
    internalTimesStoredCounter = 0;
    lastEvaluatedSpeceisForNewReactions = 0;
    nonCatalyticLastID = 0;
    // TO BE PARAMETRIZED
    Currentattempts = 0;
    resetReactionsCounter();

    if(debugLevel == FINDERRORDURINGRUNTIME) cout << "environment::environment end" << endl;

}

/* ************************
 |
 |	CREATION FUNCTIONS
 |
 *************************/

/**
 Initial molecule population creation. If the number of species stored in the configuration file is grater than the possible 
 number of species according to the alphabet and maximum length all species up to the MAX length will be created
 @version 1.0
 @param MTRand& tmpRndDoubleGen
 */
//TR void environment::createInitialMoleculesPopulation(MTRand& tmpRndDoubleGen)
//{
//	if(debugLevel == FINDERRORDURINGRUNTIME) cout << "environment::createInitialMoleculesPopulation start" << endl;
//
//	unsigned int alphabetLength = alphabet.length();
//	totalPossibleNumberOfSpecies = 0;
//	//TOTAL NUMBER OF POSSIBLE SPECIES COMPUTATION
//	int si; //Signed i
//	for(unsigned int i = 1; i <= initialMaxLength; ++i)
//	{
//		si = i;
//		totalPossibleNumberOfSpecies += pow(acs_double(alphabetLength),si);
//	}
//	//PROBABILITY LENGTH SPECIES COMPUTATION
//	const int LENGTHS = initialMaxLength;
//	acs_double proportionalInitialLengthProbability[LENGTHS];
//	for(unsigned int i = 0; i < initialMaxLength; ++i)
//	{
//		si=i;
//		if(i == 0)
//			proportionalInitialLengthProbability[i] = pow(acs_double(alphabetLength),si+1) / totalPossibleNumberOfSpecies;
//		else
//			proportionalInitialLengthProbability[i] = proportionalInitialLengthProbability[i-1] + (pow(acs_double(alphabetLength),si+1) / totalPossibleNumberOfSpecies);	
//		
//	}
//	
//	
//	//acs_int initialAmount;						// initial species amount
//	acs_double initialConcentration;			// initial concentration
//	acs_int newSpeciesLenght;					// length of the new molecule
//	acs_int monomerID;							// monomer identificator during cycles
//	acs_longInt speciesID;						// species identificator during cycles
//	acs_int fixedSpeciesID;						// fixed species identificator during the creation of monomers species
//	acs_int speciesMonomerPos;					// identificator or the monomer position inside the species during the single species creation
//	string newSpeciesSequence;					// temp variable containing the new species sequence creation
//	acs_double tmpDiffusionContribute;			// Diffusion renforcmente contribute. This value is proportional to the length of the species	
//	acs_int tmpSolubility;		// Precipitation contribute. This value is proportional to the length of the the species
//	
//	// MOLECULES CREATION
//		// IN ANY CASE MONOMERS ARE CREATED AS SPECIES
//	for(fixedSpeciesID = 0; fixedSpeciesID < alphabetLength; fixedSpeciesID++)
//	{
//		newSpeciesSequence = alphabet.substr(fixedSpeciesID, 1);
//		initialConcentration = createInitialConcentration(alphabetLength, 1);
//		// PRECIPITATION AND DIFFUSION CONTRIBUTE COMPUTATION
//		tmpDiffusionContribute = createDiffusionRenforcement(diffusion_contribute, 1);
//                tmpSolubility = setSolubility(1, tmpRndDoubleGen);
//		allSpecies.push_back(species(fixedSpeciesID, newSpeciesSequence, initialConcentration, tmpDiffusionContribute, tmpSolubility, 0, 0, 0, volume, K_nrg));
//                incNumberOfSpecies(); // incSpeciesProcedure(fixedSpeciesID) cannot be use due to the initialConcentration != 1 molecule
//		numberOfMolecules += allSpecies.at(fixedSpeciesID).getAmount();
//	}
//	
//	// ALL THE OTHERS MOLECULES WILL BE CREATED
//		// IF THE NUMBER OF INITIAL SPECIES IS LESS THAN THE NUMBER OF POSSIBILE SPECIES SPECIES AND SPECIES LENGTH WILL BE CREATED RANDOMLY
//		// OTHERWISE ALL SPECIES ACCORDING TO ALPHABET AND MAXIMUM LENGHT WILL BE CREATED
//	if(debugLevel >= RUNNING_VERSION){cout << "|- TOTAL CONCEIVABLE NUMBER OF SPECIES: " << totalPossibleNumberOfSpecies << endl;}
//	
//	if(initialPopulationNumber < totalPossibleNumberOfSpecies)
//	{
//		if(debugLevel >= RUNNING_VERSION){cout << "|- The system is going to create " << initialPopulationNumber << " species" << endl;}
//		
//		bool flagEqualSeqBetweenTwoSpecies;	// flag variable that check whether the sequence is already present
//		bool flagSpeciesStored;				// boolean variable necessary to stop while iteration
//		const long int ENDEAVOURS = 100;			// Max number of temptative to perform before to compute the next molecule if the molecule is already present
//		int endeavourCounter;				// Endeavours counter
//		
//		acs_longInt realSpeciesID = fixedSpeciesID; //THIS COUNTER IS CREATED IN RELATION TO THE POSSIBILITY TO HAVE HOLE IN NUMERATION DUE TO 
//													//THE CREATION OF ALREADY PRESENT SPECIES
//		for(speciesID = fixedSpeciesID; speciesID < initialPopulationNumber; speciesID++)
//		{
//			if(debugLevel >= SMALL_DEBUG)
//				cout << "Species Number " << speciesID << endl;
//			flagSpeciesStored = false;
//			endeavourCounter = 1;
//			flagEqualSeqBetweenTwoSpecies = true;
//			while(flagSpeciesStored == false && endeavourCounter <= ENDEAVOURS)
//			{
//				//SELECT RANDOMLY ACCORDING TO WEIGHT PROBABILITY THE LENGTH OF THE NEW MOLECULE
//				//The value of the function is incremented because the function returns the position of the weight vector
//				newSpeciesLenght = returnSelectionIdFromAWeightProbVector(proportionalInitialLengthProbability, tmpRndDoubleGen) + 1;
//				//CREATE RANDOM MOLECULE SEQUENCE
//				monomerID = getIntRandom(0, alphabetLength-1, tmpRndDoubleGen);
//				newSpeciesSequence = alphabet.substr(monomerID,1);
//				if (newSpeciesLenght > 1)
//				{
//					for(speciesMonomerPos = 2; speciesMonomerPos <= newSpeciesLenght; speciesMonomerPos++)
//					{
//						//Create sequence
//						monomerID = getIntRandom(0, alphabetLength-1, tmpRndDoubleGen);
//						newSpeciesSequence += alphabet.substr(monomerID, 1);
//					}
//				}
//				//CHECK WHETHER THE SPECIES IS ALREADY PRESENT
//				for(acs_longInt i = 0; i < (acs_longInt)allSpecies.size(); i++)	
//				{
//					if(allSpecies.at(i).getSequence() == newSpeciesSequence)
//					{
//						flagEqualSeqBetweenTwoSpecies = false;
//						if(debugLevel >= 3)
//							cout << endeavourCounter << " " << i << " " << allSpecies.at(i).getSequence() << " " << newSpeciesSequence << endl;
//						endeavourCounter++;
//						break;
//					}else{
//						flagEqualSeqBetweenTwoSpecies = true; //Otherwise the system searching for the equal species until ENDEAVOURS will be reached
//					}
//				}		
//				//CREATE NEW MOLECULE OBJECT AND STORE IT INTO THE ENVIRONMENT SPECIES VECTOR
//				if(flagEqualSeqBetweenTwoSpecies == true)
//				{
//					// The initial amount is inversaly proportional to the length
//					initialConcentration = createInitialConcentration(alphabetLength, newSpeciesLenght);					
//					// PRECIPITATION AND DIFFUSION CONTRIBUTE COMPUTATION
//					tmpDiffusionContribute = createDiffusionRenforcement(diffusion_contribute, newSpeciesLenght);
//                                        tmpSolubility = setSolubility(newSpeciesLenght, tmpRndDoubleGen);
//					allSpecies.push_back(species(realSpeciesID, newSpeciesSequence, initialConcentration, tmpDiffusionContribute, tmpSolubility, 0, 0, 0, volume, K_nrg));
//					flagSpeciesStored = true;
//                                        if(initialConcentration > 0)
//                                            numberOfSpecies++;
//					numberOfMolecules += allSpecies.at(realSpeciesID).getAmount();
//					realSpeciesID++;
//				}
//			}
//		}
//		
//		if(debugLevel >= RUNNING_VERSION){cout << "\t|-" << (acs_longInt)allSpecies.size() << " SPECIES HAS BEEN CREATED" << endl;}
//		
//	}else{
//		// ALL POSSIBILE SPECIES CREATION
//		if(debugLevel == RUNNING_VERSION){ cout << "|- ALL POSSIBILE SPECIES CREATION (initial number of species has been set greater than the conceivable number)" << endl;}
//		
//		//initialize temporary Vector to store species population
//		vector<species> tmpAllSpecies;
//		// UPdate species ID considering monomers
//		speciesID = fixedSpeciesID;
//		// Declare boolean remove flag variable and removedSpeciesLengths iterator
//		int removedSpeciesLengths;
//		for(unsigned int tmpSpeciesLength = 2; tmpSpeciesLength <= initialMaxLength; tmpSpeciesLength++)
//		{ //For each possibile initial species length
//			//CREATE A TEMPORARY COPY OF THE SPECIES VECTOR
//			tmpAllSpecies = allSpecies;
//			//REMOVE ALL NOT NECESSARY (THE ONES WITH LENGTH < tmpSpeciesLength - 1) SPECIES FROM THE TEMPORARY VECTOR
//			removedSpeciesLengths = 1;
//			if (tmpSpeciesLength > 1)
//			{
//				while (tmpAllSpecies.at(0).getSequence().length() < (tmpSpeciesLength - 1))
//				{
//					    tmpAllSpecies.removeFirst();
//				}
//			}
//			
//			for(unsigned int tmpAlphabetPosition = 0; tmpAlphabetPosition < alphabetLength; tmpAlphabetPosition++)
//			{ //for each possibile monomer
//				if(debugLevel >= HIGH_DEBUG){cout << "\t\t\t|- Monomer to add -> " << alphabet[tmpAlphabetPosition] << endl;}
//				
//				// TO CRATE ALL POSSIBILE SPECIES UP TO LENGTH M RECURSIVELY MONOMERS WILL BE ADDED AT THE END OF THE PRESENT SPECIES
//				for(acs_longInt i = 0; i < (acs_longInt)tmpAllSpecies.size(); i++)
//				{
//					newSpeciesSequence = tmpAllSpecies.at(i).getSequence() + alphabet[tmpAlphabetPosition];
//					if(debugLevel >= HIGH_DEBUG){cout << "\t\t\t\t" << tmpAllSpecies.at(i).getSequence() << 
//														" + " << alphabet[tmpAlphabetPosition] << " = " << 
//														newSpeciesSequence << endl;}
//					
//					// The initial amount is inversaly proportional to the length
//					newSpeciesLenght = newSpeciesSequence.length();
//					initialConcentration = createInitialConcentration(alphabetLength, newSpeciesLenght);	
//					// PRECIPITATION AND DIFFUSION CONTRIBUTE COMPUTATION
//					tmpDiffusionContribute = createDiffusionRenforcement(diffusion_contribute, newSpeciesLenght);
//                                        tmpSolubility = setSolubility(newSpeciesLenght, tmpRndDoubleGen);
//					allSpecies.push_back(species(speciesID, newSpeciesSequence, initialConcentration, tmpDiffusionContribute, 
//												 tmpSolubility, 0, 0, 0, volume, K_nrg));
//                                        if(initialConcentration > 0)
//                                            numberOfSpecies++;
//					numberOfMolecules += allSpecies.at(speciesID).getAmount();
//					speciesID++; //update species ID
//				}				
//			}
//		}
//	}
//	
//	// FIRING DISK STORING, INFLUX CREATION AND PROBABILITY VECTOR FOR A NUTRIENT TO BE SELECTED IS CREATED
//	firingDisk = allSpecies;
//	this->getNutrientsFromTheFiringDisk();
//	if(debugLevel >= SMALL_DEBUG){printNutrientsAndProbability();}
//	if((initialMaxLength > influx) & (influx > 0))
//	{
//		this->nutrientsAmountsFixing();
//	}
//	
//	if(debugLevel >= SMALL_DEBUG){printNutrientsAndProbability();}
//	
//	if(debugLevel >= RUNNING_VERSION){cout << "\t\t\t|- DONE!!! " << endl;}
//	
//	if(debugLevel == FINDERRORDURINGRUNTIME) cout << "environment::createInitialMoleculesPopulation end" << endl;
//
//}// eof createInitialMoleculesPopulation

/**
 initial layer initialization
 @version 1.0
 @param MTRand& tmpRndDoubleGen
*/
//void environment::createReactionsLayer(MTRand& tmpRndDoubleGen)
//{
//	if(debugLevel == FINDERRORDURINGRUNTIME) cout << "environment::createReactionsLayer start" << endl;
//
//	if(debugLevel >= RUNNING_VERSION){cout << "|- Reactions Layer creation..." << endl;}
//	
//	// COPY ALL THE SPECIES IN A TEMPORARY VECTOR
//	vector<acs_longInt> allPresentSpecies;
//
//	for(acs_longInt i = 0; i < (acs_longInt)allSpecies.size(); i++)
//	{
//			allPresentSpecies.push_back(i);
//	}
//	
//	// COMPUTE TOTAL NUMBER OF CONCEIVABLE REACTIONS
//	acs_longInt totalNumberOfConcevableReactions = 0;
//	if((acs_longInt)allPresentSpecies.size() > 0)
//	{
//		for(acs_longInt sID = 0; sID < (acs_longInt)allPresentSpecies.size(); sID++)
//		{
//			totalNumberOfConcevableReactions += allSpecies.at(allPresentSpecies.at(sID)).getSequenceLength() - 1;
//		}	
//	
//		totalNumberOfConcevableReactions += pow((double)allPresentSpecies.size(), 2.0);
//	}
//	
//	if(debugLevel >= SMALL_DEBUG){printAllSpeciesIdAndSequence();}
//	
//	if(debugLevel >= RUNNING_VERSION){cout << "|- Total number of conceivable reactions: " << totalNumberOfConcevableReactions << endl;}
//	
//	//FOR EACH SPECIES
//	acs_int reactionsForThisSpecies;
//	for(acs_longInt sID = 0; sID < (acs_longInt)allPresentSpecies.size(); sID++)
//	{
//		
//		// Only species longer than nonCatalyticMaxLength can catalyse reactions!!!
//		if(allSpecies.at(allPresentSpecies.at(sID)).getSequenceLength() > nonCatalyticMaxLength)
//		{
//			reactionsForThisSpecies = computeSngSpeciesRctsNumber(totalNumberOfConcevableReactions, tmpRndDoubleGen);
//		}else{
//			reactionsForThisSpecies = 0;
//		}
//		
//		// CREATE REACTIONS FOR THIS SPECIFIC SPECIES
//		if(debugLevel >= RUNNING_VERSION){cout << "\t|- Creating reactions for species : " << sID << "[" << allSpecies.at(sID).getSequence() << "]" << endl;}
//		if(reactionsForThisSpecies > 0)
//		{
//			if(!createReactionsForThisSpecies(sID, reactionsForThisSpecies, tmpRndDoubleGen, allPresentSpecies, NEWREACTIONS))
//			{	
//				ExitWithError("createReactionsForThisSpecies", "	|- !*!*!*! Problem with the reactions creation");
//			}else{
//				allSpecies[sID].setEvaluated(); // IF the reactions have been created the molecules evalutation is setting on 1
//			}
//		}else{
//			allSpecies[sID].setEvaluated(); // IF the reactions have been created the molecules evalutation is setting on 1		
//		}
//	}
//	if(debugLevel >= RUNNING_VERSION){cout << "\t\t\t|- DONE!!! " << endl;}
//
//	if(debugLevel == FINDERRORDURINGRUNTIME) cout << "environment::createReactionsLayer end" << endl;
//
//}//eof

/**
 Compute number of reaction catalysd by a catalyst according to the total number of reactions and reactions probabilities
 @version 1.0
 @param acs_int tmpTotalNumberOfReactions Total number of conceivable reactions
 @param acs_double tmpRctsProb reaction probability
 */
acs_int environment::computeSngSpeciesRctsNumber(acs_longInt tmpTotalNumberOfReactions, MTRand& tmpRndDoubleGen)
{
	if(debugLevel == FINDERRORDURINGRUNTIME) cout << "\t\t\t\tenvironment::computeSngSpeciesRctsNumber start" << endl;

	acs_int temp_tmpReactionsForThisSpecies; 
	acs_int tmpReactionsForThisSpecies; 
	
	temp_tmpReactionsForThisSpecies = tmpTotalNumberOfReactions * reactionProbability;
	
	if(debugLevel >= SMALL_DEBUG)
		cout << "!MESSAGE :: Reaction for this species BEFORE stochastic fluctuation: " << temp_tmpReactionsForThisSpecies << endl;
	
	if(temp_tmpReactionsForThisSpecies <= 100)
	{
		//tmpReactionsForThisSpecies = random_poisson((acs_double)temp_tmpReactionsForThisSpecies	, tmpRndDoubleGen);
		tmpReactionsForThisSpecies = random_binomial(tmpTotalNumberOfReactions, reactionProbability, tmpRndDoubleGen);
	}else {
		tmpReactionsForThisSpecies = temp_tmpReactionsForThisSpecies;
	}
		
	if(tmpReactionsForThisSpecies <= 0)
		tmpReactionsForThisSpecies = 0;
	
	if(debugLevel >= SMALL_DEBUG)
		cout << "!MESSAGE :: Reaction for this species AFTER stochastic fluctuation: " << tmpReactionsForThisSpecies << endl;
		
	if(debugLevel == FINDERRORDURINGRUNTIME) cout << "\t\t\t\tenvironment::computeSngSpeciesRctsNumber end" << endl;

	return tmpReactionsForThisSpecies;
	

}//eof


/**
 Creation of all the reactions related to one specific species
 @version 1.1
 @date 2011/07/07
 @param acs_longInt tmpsID species vector ID
 @param acs_int tmpReactionsForThisSpecies number of reactions to create for this species
 @param MTRand& tmp_RndDoubleGen random number generator
 @param vector<acs_longInt>& tmpIDOfCandidateSpecies ID of the species avalaible for the reaction
 @param acs_int tmpRctCreationType NEWREACTION or UPGRADEREACTIONS
 */
bool environment::createReactionsForThisSpecies(acs_longInt tmpsID, acs_int tmpReactionsForThisSpecies, MTRand& tmp_RndDoubleGen, 
                                                vector<acs_longInt>& tmpIDOfCandidateSpecies, acs_int tmpRctCreationType)
{	
    if(debugLevel == FINDERRORDURINGRUNTIME) cout << "\t\t\t\tenvironment::createReactionsForThisSpecies start" << endl;

    acs_int reactionType;
    acs_longInt ids_I = 0;
    acs_longInt ids_II = 0;
    acs_longInt ids_III = 0;
    acs_double tmpK_ass; // Declare three kinetic variable to initialize according to the presence of not of the reverse reations
    acs_double tmpK_diss;
    acs_double tmpK_cpx;
    acs_double tmpK_spontanoues; // This variable is used to set the spontanoues kinetic constant of the reaction (0 in case of no spontaneous)
    acs_longInt id_reaction;
    acs_int tempEnergyRctType = ESOERGONIC;
    acs_int nrgRctBoolFncID;
    acs_int tmpEnergizable = NOTENERGIZABLE;
    bool validReactionFlag = false; // this flag is true if all possible steps of the reaction creation have been satisfied
    bool toProceed = true; // Variable to check whether it is possibile proceeding with the creation reaction or not
    acs_int rctCreationAttemptCounter = 0; // Number of tentative to create a reaction
	
    // CREATE A TEMPORARY VECTOR CONTAINING ALL USEFUL SPECIES
    vector<acs_longInt> speciesAvailableForReactions;
    vector<acs_longInt> speciesAvailableForCleavageReactions;
	
    for(acs_longInt i = 0; i < (acs_longInt)tmpIDOfCandidateSpecies.size(); i++)
    {
       try{
            if(allSpecies.at(tmpIDOfCandidateSpecies.at(i)).getComplexCutPnt() == 0)
            {
				if(tmpRctCreationType == NEWREACTIONS) // If this is a reactions for a new species
				{
					speciesAvailableForReactions.push_back(tmpIDOfCandidateSpecies.at(i));
				}else{
					// OLD SPECIES REACTIONS UPDATING PROCESS CAN INVOLVE SPECIES NEVER CONSIDERED BEFORE BY THIS SPECIES
					// Species belonging to the firing disk are not considered
					if((tmpIDOfCandidateSpecies.at(i) > ((acs_int)firingDisk.size() - 1)) &&
					   (tmpIDOfCandidateSpecies.at(i) > tmpsID) &&
					   (tmpIDOfCandidateSpecies.at(i) > allSpecies.at(tmpsID).getLastSpeciesEvaluated()))
							speciesAvailableForReactions.push_back(tmpIDOfCandidateSpecies.at(i));
				}
            }else{
                    ExitWithError("createReactionsForThisSpecies", "A Complex Cannot be envolved in a reaction as product or substrate");
            }
        }catch(exception&e)
        {
            cout << "createReactionsForThisSpecies allSpecies.at(tmpIDOfCandidateSpecies.at(i)).getComplexCutPnt()" << endl;
            cout << "Vectorsize " << tmpIDOfCandidateSpecies.size()<<" - position " << i << endl;
            cerr << "exceptioncaught:" << e.what() << endl;
            ExitWithError("createReactionsForThisSpecies","exceptionerrorthrown");
        }

    }

    if((acs_longInt)speciesAvailableForReactions.size() > 0) // If at least one species may be interact with this species
    {
    // COPY AVALAIBLE SPECIES VECTOR INTO AVALAIBLE FOR CLEAVAGE
        speciesAvailableForCleavageReactions = speciesAvailableForReactions;
        allSpecies.at(tmpsID).setLastSpeciesEvaluated(speciesAvailableForReactions.at(speciesAvailableForReactions.size()-1));
        if(debugLevel >= SMALL_DEBUG)
        {
                cout << "\t\t|- " << allSpecies.at(tmpsID).getID() + 1 << "th Species (ID " << allSpecies.at(tmpsID).getID() << ") - " <<
                (acs_longInt)tmpIDOfCandidateSpecies.size() << " Candidate substrates - Number of Reactions: " << tmpReactionsForThisSpecies <<
                " - Total Reactions: " << (acs_longInt)allReactions.size() <<  endl;
        }

        if(debugLevel >= SMALL_DEBUG)
        {
                cout << "\t ----------------------------------------------" << endl;
                cout << "\t ----------------------------------------------" << endl;
        }
        for(acs_int IDrct = 0; IDrct < tmpReactionsForThisSpecies; IDrct++) // FOR EACH REACTION
        {
            if(debugLevel >= SMALL_DEBUG)
            {
                    cout << "\t ----------------------------------------------" << endl;
                    cout << "\t|- Species " << tmpsID << " - Reaction Number " << IDrct <<
                    " of " << tmpReactionsForThisSpecies << endl;
            }
            rctCreationAttemptCounter = 0;
            validReactionFlag = false;
            while(!validReactionFlag && (rctCreationAttemptCounter < 200))
            {
                    tempEnergyRctType = ESOERGONIC;
                    validReactionFlag = false; //RESET VALIDATION FLAG IN ORDER TO NOT CREATE WRONG REACTIONS
                    toProceed = true;

                    if(rctCreationAttemptCounter == 0)
                    {
                            // If this is the first temptative to create a reaction reaction direction is randomly choosen
                            reactionType = selectWhetherCleavageOrCond(tmp_RndDoubleGen);
                    }

                    //------------------------------ CLEAVAGE REACTION -------------------------------------------------|

                    if(reactionType == CLEAVAGE)
                    {
                       acs_longInt substratePosID;
                       if((acs_longInt)speciesAvailableForCleavageReactions.size() > 0) //IF THERE ARE SPECIES AVALAIBLE for CLEAVAGE
                        {
                            // Select species to hydrolize
                           substratePosID = returnUniformSelection_LONG_IdFromVector(speciesAvailableForCleavageReactions, tmp_RndDoubleGen);
                           acs_int substrateLength;
                           try{
                               substrateLength = allSpecies.at(speciesAvailableForCleavageReactions.at(substratePosID)).getSequence().length();
                           }catch(exception&e){
                               cout << "createReactionsForThisSpecies allSpecies.at(speciesAvailableForCleavageReactions.at(substratePosID))" << endl;
                               cout << "Vectorsize "<<allSpecies.size()<<" - position "<<speciesAvailableForCleavageReactions.at(substratePosID)<<endl;
                               cerr << "exceptioncaught:"<<e.what()<<endl;
                               ExitWithError("createReactionsForThisSpecies","exceptionerrorthrown");
                           }

                           string product_I;
                           string product_II;
                           if(substrateLength > 1) //IF THE SUBSTRATE IS LONGER THAN 1
                           {
                               if(debugLevel >= SMALL_DEBUG)
                               {
                                    //cout << "\t\t -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-" << endl;
                                    //cout << "\t\t|- Reaction: " << IDrct << " CLEAVAGE REACTION:" << endl;
                                    //cout << "\t\t\t|- SUBSTRATE pos: " << substratePosID << " | SUBSTRATE ID: " <<
                                    //allSpecies.at(speciesAvailableForCleavageReactions.at(substratePosID)).getID() << " | SEQUENCE: " <<
                                    //allSpecies.at(speciesAvailableForCleavageReactions.at(substratePosID)).getSequence() << endl;
                                    //cout << "\t\t\t|- Total number of avalaible cleavage: " << (acs_longInt)speciesAvailableForCleavageReactions.size() << endl;
                               }
                                //SELECT CUTTING POINT AND CREATE THE TWO PRODUCTS
                                acs_int cuttingPoint = getIntRandom(1, substrateLength - 1, tmp_RndDoubleGen);

                                try{
                                    product_I = allSpecies.at(speciesAvailableForCleavageReactions.at(substratePosID)).getSequence().substr(0, cuttingPoint);
                                    product_II = allSpecies.at(speciesAvailableForCleavageReactions.at(substratePosID)).getSequence().substr(cuttingPoint, allSpecies.at(speciesAvailableForCleavageReactions.at(substratePosID)).getSequence().length() - cuttingPoint);
                                }catch(exception&e){
                                    cout << "createReactionsForThisSpecies allSpecies.at(speciesAvailableForCleavageReactions.at(substratePosID))" << endl;
                                    cout << "Vectorsize "<<allSpecies.size()<<" - position "<<speciesAvailableForCleavageReactions.at(substratePosID)<<endl;
                                    cerr << "exceptioncaught in if(substrateLength > 1):"<<e.what()<<endl;
                                    ExitWithError("createReactionsForThisSpecies","exceptionerrorthrown");
                                }

                                if(debugLevel >= MEDIUM_DEBUG)
                                {
                                    cout << "\t\t\t\t|- Product 1: " << product_I << endl;
                                    cout << "\t\t\t\t|- Product 2: " << product_II << endl;
                                }

                                // Try here is not necessary since it has been used 20 rows above
                                ids_I = speciesAvailableForCleavageReactions.at(substratePosID); //substrate vector allSpecies position

                                // CHECK WHETHER THE NEW PRODUCT IS ALREADY PRESENT AND IF IT ISN'T REACTION IS STORED. IF IT IS, ITS ID BECAME THE VALUE TO STORE INTO THE REACTIONS,
                                // OTHERWISE A NEW ID WILL BE CREATED
                                ids_II = returnPosSpeciesAlreadyPresent(product_I);
                                if(product_I != product_II) // IF PRODUCTS ARE DIFFERENT
                                {
                                    ids_III = returnPosSpeciesAlreadyPresent(product_II);
                                    // IF BOTH PRODUCTS ARE NEW THE ID OF HTE SECOND ONE IN INCREMENTED BY 1
                                    if(ids_II == ids_III)
                                            ids_III++;
                                }else{
                                    ids_III = ids_II;
                                }

                                // This code segment has been commented since the introduction of the external initializator.
                                // SELF REPLICATION IS AVOIDED INTO THE FIRING DISK, IF ONE OF THE PRODUCTS IS A CATALYST (IN THE FIRING DISK) THE REACTION IS NOT POSSIBLE
//                                  if(tmpsID <= ((acs_int)firingDisk.size() - 1))
//                                  {
//                                      if(tmpsID == ids_II) toProceed = false;
//					if(tmpsID == ids_III) toProceed = false;
//
//					// CHECK IF THE CANDIDATE PRODUCTS ALREADY CATALYSE THE CATALYST IN ORDER TO AVOID MUTUAL CATALYSIS
//					if(toProceed) toProceed = checkIfOnlyMutualCatalysis(tmpsID, ids_II);
//                                      if(toProceed && (ids_II != ids_III)) toProceed = checkIfOnlyMutualCatalysis(tmpsID, ids_III);
//                                  }

                                if(debugLevel >= MEDIUM_DEBUG)
                                {
                                        cout << "\t\t\t\t|- ID 1: " << ids_II << endl;
                                        cout << "\t\t\t\t|- ID 2: " << ids_III << endl;
                                }

                                // ONE SPECIES CAN NOT CATALYZE BOTH FORWARD AND BACKWARD REACTIONS
                                if(toProceed)
                                {
                                    if(!notInverseReactionAlreadyCatalyzed(CLEAVAGE, ids_I, ids_II))
                                    {
                                        validReactionFlag = false; // BACKWARD REACTION IS ALREADY CATALIZED!!!
                                        if(debugLevel == SMALL_DEBUG)
                                                cout << "\t\t\tFound COND reverse reaction - " << validReactionFlag <<  " -  tentativi - " << rctCreationAttemptCounter <<
                                                " - Numero " << IDrct << endl;
                                    }else{
                                      if(ids_II == (acs_longInt)allSpecies.size()) //IF THE POSITION RETURNED IS EQUAL TO THE SIZE OF vector ALLSPECIES, IT MEANS THAT THE SPECIES IS NOT PRESENT
                                        {
                                                if(tmp_RndDoubleGen() < ratioSpeciesEnergizable){tmpEnergizable=ENERGIZABLE;}else{tmpEnergizable=NOTENERGIZABLE;}
                                                allSpecies.push_back(species(ids_II, product_I, (acs_longInt)0, 0, 0, 0, 0, 0, volume, 0, tmpEnergizable, influx_rate, maxLOut));
                                        }
                                      if(ids_III == (acs_longInt)allSpecies.size()) //IF THE POSITION RETURNED IS EQUAL TO THE SIZE OF vector ALLSPECIES, IT MEANS THAT THE SPECIES IS NOT PRESENT
                                        {
                                                if(tmp_RndDoubleGen() < ratioSpeciesEnergizable){tmpEnergizable=ENERGIZABLE;}else{tmpEnergizable=NOTENERGIZABLE;}
                                                allSpecies.push_back(species(ids_III, product_II, (acs_longInt)0, 0, 0, 0, 0, 0, volume, 0, tmpEnergizable, influx_rate, maxLOut));
                                        }
                                        validReactionFlag = true;
                                    }
                                }else{
                                        validReactionFlag = false;
                                }
                            }
                            // REMOVE THIS SUBSTRATE FOR THE AVALAIBLE SUBSTRATE
                            // speciesAvailableForCleavageReactions.removeAt(substratePosID);
                            speciesAvailableForCleavageReactions.erase(speciesAvailableForCleavageReactions.begin()+substratePosID);
                        }

                    //------------------------------ CONDENSATION REACTION -------------------------------------------------|

                    }else{
                        acs_longInt substratePosID_I;
                        acs_longInt substratePosID_II;
                        string product;
                        if((acs_longInt)speciesAvailableForReactions.size() > 0) // IF THERE ARE SPECIES
                        {
                            substratePosID_I = returnUniformSelection_LONG_IdFromVector(speciesAvailableForReactions, tmp_RndDoubleGen);
                            substratePosID_II = returnUniformSelection_LONG_IdFromVector(speciesAvailableForReactions, tmp_RndDoubleGen);

                            if(debugLevel >= SMALL_DEBUG)
                            {
                                cout << "\t\t -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-" << endl;
                                cout << "\t\t|- Reaction: " << IDrct << " CONDENSATION REACTION:" << endl;
                                cout << "\t\t\t|- SUBSTRATE 1 ID: " << allSpecies.at(speciesAvailableForReactions.at(substratePosID_I)).getID()
                                         << " | SEQUENCE: " << allSpecies.at(speciesAvailableForReactions.at(substratePosID_I)).getSequence() << endl;
                                cout << "\t\t\t|- SUBSTRATE 2 ID: " << allSpecies.at(speciesAvailableForReactions.at(substratePosID_II)).getID()
                                         << " | SEQUENCE: " << allSpecies.at(speciesAvailableForReactions.at(substratePosID_II)).getSequence() << endl;
                            }

                            //PRODUCT CREATION
                            try{
                                product = allSpecies.at(speciesAvailableForReactions.at(substratePosID_I)).getSequence() +
                                                        allSpecies.at(speciesAvailableForReactions.at(substratePosID_II)).getSequence();
                            }catch(exception&e){
                                cout << "allSpecies.at(speciesAvailableForReactions.at(substratePosID_I)).getSequence() + allSpecies.at(speciesAvailableForReactions.at(substratePosID_II)).getSequence();" << endl;
                                cout << "Vectorsize "<<allSpecies.size()<<" - position "<<speciesAvailableForCleavageReactions.at(substratePosID_I)<<endl;
                                cerr << "exceptioncaught in if(substrateLength > 1):"<<e.what()<<endl;
                                ExitWithError("createReactionsForThisSpecies","exceptionerrorthrown");
                            }

                            if(debugLevel >= SMALL_DEBUG)
                                cout << "\t\t\t\t|- Product: " << product << endl;

                            ids_I = returnPosSpeciesAlreadyPresent(product);

                            if(tmpsID <= ((acs_int)firingDisk.size() - 1))
                            {
                                if(tmpsID == ids_I) toProceed = false;

                                // CHECK WHETHER THE CANDIDATE PRODUCT ALREADY CATALYSES THE CATALYST
//							if(toProceed) toProceed = checkIfOnlyMutualCatalysis(tmpsID, ids_I);
                            }

                            if(toProceed)
                            {
                                // ONE SPECIES CAN NOT CATALYZE BOTH ONE REACTION AND THE INVERTED ONE
                                try{
                                    ids_II = speciesAvailableForReactions.at(substratePosID_I);
                                }catch(exception&e){
                                    cout << "speciesAvailableForReactions.at(substratePosID_I)" << endl;
                                    cout << "Vectorsize "<<speciesAvailableForReactions.size()<<" - position "<<substratePosID_I<<endl;
                                    cerr << "ids_II = speciesAvailableForReactions.at(substratePosID_I);"<<e.what()<<endl;
                                    ExitWithError("createReactionsForThisSpecies","exceptionerrorthrown");
                                }

                                try{
                                    ids_III = speciesAvailableForReactions.at(substratePosID_II);
                                }catch(exception&e){
                                    cout << "speciesAvailableForReactions.at(substratePosID_II)" << endl;
                                    cout << "Vectorsize "<<speciesAvailableForReactions.size()<<" - position "<<substratePosID_II<<endl;
                                    cerr << "ids_II = speciesAvailableForReactions.at(substratePosID_II);"<<e.what()<<endl;
                                    ExitWithError("createReactionsForThisSpecies","exceptionerrorthrown");
                                }

                                if(!notInverseReactionAlreadyCatalyzed(CONDENSATION, ids_I, ids_II))
                                {
                                    validReactionFlag = false;
                                    if(debugLevel == SMALL_DEBUG)
                                         cout << "\t\t\tFound clavage reverse reaction - " << validReactionFlag <<  " -  tentativi - " << rctCreationAttemptCounter <<
                                         " - Numero " << IDrct << endl;
                                }else{
                                    if(ids_I == (acs_longInt)allSpecies.size()) //IF THE POSITION RETURNED IS EQUAL TO THE SIZE OF VECOTR ALLSPECIES, IT MEANS THAT THE SPECIES IS NOT PRESENT
                                    {
                                            if(tmp_RndDoubleGen() < ratioSpeciesEnergizable){tmpEnergizable=ENERGIZABLE;}else{tmpEnergizable=NOTENERGIZABLE;}
                                            allSpecies.push_back(species(ids_I, product, (acs_longInt)0, 0, 0, 0, 0, 0, volume, 0, tmpEnergizable, influx_rate, maxLOut));
                                    }
                                    validReactionFlag = true;
                                }
                            }else{
                                validReactionFlag = false;
                            }
                        }

                    }

                    if(debugLevel == SMALL_DEBUG)
                        cout << "\t\t\tReaction validity: " << validReactionFlag << " - attempts " << rctCreationAttemptCounter << endl;

                    if(validReactionFlag == true)
                    {
                        // IF ALL IT'S OK REACTION AND CATALYSIS ARE STORED
                        acs_int tmpCpxTarget = 1;
    					if(tmp_RndDoubleGen() >= 0.5) tmpCpxTarget = 2;

                        // SINCE THE SPECIES HAS BEEN CREATED, NOW I HAVE TO CHECK WHETHER THE REACTION IS ALREADY PRESENT.
                        if(getNumberOfReactions() > 0)
                        {
                            id_reaction = returnPosReactionAlreadyPresent(reactionType, ids_I, ids_II, ids_III);
                        }else{
                            id_reaction = 0;
                        }

                        if(id_reaction == getNumberOfReactions())
                        {
                            if(debugLevel == SMALL_DEBUG)
                                cout << "\t\t\tStoring Reactions and Assagning catalysis for reaction " << IDrct << " - ID REACTION STORED: " << id_reaction << endl;
                            // IF THE ENERGY IS CONSIDERED
                            // According to the energy parameter and the reaction type the energetic type is assigned
                            if(nrgBoolFlag == ENERGYBASED)
                            {
                                if(nrgBooleanFunctions.size() > 0)
                                {
                                    nrgRctBoolFncID = returnSelectionIdFromAWeightProbVectorAlreadyNormalized(nrgBoolFncsProb2BeSelected,tmp_RndDoubleGen);
                                    try{
                                        tempEnergyRctType = nrgBooleanFunctions.at(nrgRctBoolFncID);
                                    }catch(exception&e){
                                        cout << "tempEnergyRctType = nrgBooleanFunctions.at(nrgRctBoolFncID);" << endl;
                                        cout << "Vectorsize "<<nrgBooleanFunctions.size()<<" - position "<<nrgRctBoolFncID<<endl;
                                        cerr << "tempEnergyRctType = nrgBooleanFunctions.at(nrgRctBoolFncID);"<<e.what()<<endl;
                                        ExitWithError("createReactionsForThisSpecies","exceptionerrorthrown");
                                    }
                                }else{
                                    // If the vector contains only 1 element, the first one is selected
                                    tempEnergyRctType = nrgBooleanFunctions.front();
                                }
                            }
                            // IF IT IS NOT PRESENT THE VECTOR ALLREACTIONS AND THE VECTOR ALLCATALYSIS ARE UPDATED
                            if(spontRct == NOSPONTANEOUS){tmpK_spontanoues = NOSPONTANEOUS;}
                            else{
                            	if((reactionType == CONDENSATION) || (reactionType == ENDO_CONDENSATION)){tmpK_spontanoues = K_spont_ass;}
                            	else{tmpK_spontanoues = K_spont_diss;}
                            }
                            allReactions.push_back(reactions(id_reaction, reactionType, ids_I, ids_II, ids_III, 0, tempEnergyRctType, tmpK_spontanoues));
                            if(reverseReactions == true)
                            {
                                if(reactionType == CONDENSATION)
                                {
                                    tmpK_ass = K_ass;
                                    tmpK_diss = K_diss / revRctRatio;
                                    tmpK_cpx = K_cpx;
                                }else{ //Cleavage reaction
                                    tmpK_ass = K_ass / revRctRatio;
                                    tmpK_diss = K_diss;
                                    tmpK_cpx = K_cpx / revRctRatio;
                                }
                            }else{ // reverseReactions == FALSE
                                tmpK_ass = K_ass;
                                tmpK_diss = K_diss;
                                tmpK_cpx = K_cpx;
                            }

                            if(allSpecies.at(tmpsID).getSequenceLength() <= nonCatalyticMaxLength) ExitWithError("createReactionsForThisSpecies","ERROR: The catalyst is smaller than the minumum");
                            if((tmpCpxTarget < 1) || (tmpCpxTarget > 2)) ExitWithError("createReactionsForThisSpecies:r878","WRONG COMPLEX TARGET");
                            allCatalysis.push_back(catalysis(getNumberOfCatalysis(), tmpsID, id_reaction, 0, tmpK_ass, tmpK_diss, tmpK_cpx, tmpCpxTarget));

                        }else{ // OTHERWISE ONLY THE VECTOR ALLCATALYSIS IS UPDATED
                            if(reverseReactions == true)
                            {
                                if(reactionType == CONDENSATION)
                                {
                                    tmpK_ass = K_ass;
                                    tmpK_diss = K_diss / revRctRatio;
                                    tmpK_cpx = K_cpx;
                                }else{ //Cleavage reaction
                                    tmpK_ass = K_ass / revRctRatio;
                                    tmpK_diss = K_diss;
                                    tmpK_cpx = K_cpx / revRctRatio;
                                }
                            }else{ // reverseReactions == FALSE
                                tmpK_ass = K_ass;
                                tmpK_diss = K_diss;
                                tmpK_cpx = K_cpx;
                            }
                            if(debugLevel == SMALL_DEBUG)
                                cout << "\t\t\tAssagning catalysis for reaction " << IDrct << endl;
                            //Check if this reaction is not already catalysed by this species
                            if(!checkIfTheReactionIsAlreadyCatalyzedByThisSpecies(tmpsID, id_reaction))
                            {
                            	if(allSpecies.at(tmpsID).getSequenceLength() <= nonCatalyticMaxLength) ExitWithError("createReactionsForThisSpecies","ERROR: The catalyst is smaller than the minumum");
                            	if((tmpCpxTarget < 1) || (tmpCpxTarget > 2)) ExitWithError("createReactionsForThisSpecies:r905","WRONG COMPLEX TARGET");
                            	allCatalysis.push_back(catalysis(getNumberOfCatalysis(), tmpsID, id_reaction, 0, tmpK_ass, tmpK_diss, tmpK_cpx, tmpCpxTarget));
                            }
                        }
                    }else{
                            rctCreationAttemptCounter++;
                    }// if(validReactionFlag)
            } //while(!validReactionFlag)
        } // end for(acs_int IDrct = 0; IDrct < tmpReactionsForThisSpecies; IDrct++)
    }

    if(debugLevel == FINDERRORDURINGRUNTIME) cout << "\t\t\t\tenvironment::createReactionsForThisSpecies end" << endl;


    return true;
}

/**
 Old Species Reactions and Catalysis update process
 @version 1.0
 @date 2013/07/30
 @param acs_longInt tmpIDtoUpdate ID of the species to update
 @param acs_longInt ID of the new species involved in the reactions to update
 @param acs_int tmpRctType Reaction type, namely cleavage or condensation
 @param MTRand& tmp____RndDoubleGen Random number generator
 */
bool environment::updateReactions(acs_longInt tmpIDtoUpdate, acs_longInt tmpNewSpecies, acs_int tmpRctType, vector<acs_longInt>& tmp_AlreadyEvaluatedSpeciesVector, MTRand& tmp_RndDoubleGen)
{
	 if(debugLevel == FINDERRORDURINGRUNTIME) cout << "\t\t\t\tenvironment::updateReactions start" << endl;

	    acs_longInt ids_I = 0;
	    acs_longInt ids_II = 0;
	    acs_longInt ids_III = 0;
	    acs_double tmpK_ass; // Declare three kinetic variable to initialize according to the presence of not of the reverse reations
	    acs_double tmpK_diss;
	    acs_double tmpK_cpx;
	    acs_longInt id_reaction;
	    acs_double tmpK_spontanoues;
	    acs_int tempEnergyRctType = ESOERGONIC;
	    acs_int nrgRctBoolFncID;
	    acs_int tmpEnergizable = NOTENERGIZABLE;
	    bool validReactionFlag = false; // this flag is true if all possible steps of the reaction creation have been satisfied
	    bool toProceed = true; // Variable to check whether it is possibile proceeding with the creation reaction or not
	    acs_int rctCreationAttemptCounter = 0; // Number of tentative to create a reaction
		acs_int tmpCpxTarget; // Substrate involved in the complex formation (condensation reactions or reverse cleavage reactions)

	    while(!validReactionFlag && (rctCreationAttemptCounter < 200))
		{
				tempEnergyRctType = ESOERGONIC;
				validReactionFlag = false; //RESET VALIDATION FLAG IN ORDER TO NOT CREATE WRONG REACTIONS
				toProceed = true;

				//------------------------------ CLEAVAGE REACTION -------------------------------------------------|

				if(tmpRctType == CLEAVAGE)
				{
				   string product_I;
				   string product_II;
				   // Define complex target in case of reverse reaction...
				   if(tmp_RndDoubleGen() < 0.5) tmpCpxTarget = 1; else tmpCpxTarget = 2;

				   if(allSpecies.at(tmpNewSpecies).getSequenceLength() > 1) //IF THE SUBSTRATE IS LONGER THAN 1
				   {
						//SELECT CUTTING POINT AND CREATE THE TWO PRODUCTS
						acs_int cuttingPoint = getIntRandom(1, allSpecies.at(tmpNewSpecies).getSequenceLength() - 1, tmp_RndDoubleGen);

						try{
							product_I = allSpecies.at(tmpNewSpecies).getSequence().substr(0, cuttingPoint);
							product_II = allSpecies.at(tmpNewSpecies).getSequence().substr(cuttingPoint, allSpecies.at(tmpNewSpecies).getSequence().length() - cuttingPoint);
						}catch(exception&e){
							cout << "updateReactions allSpecies.at(speciesAvailableForCleavageReactions.at(substratePosID))" << endl;
							cerr << "exceptioncaught in if(allSpecies.at(tmpNewSpecies).getSequenceLength() > 1):"<<e.what()<<endl;
							ExitWithError("createReactionsForThisSpecies","exceptionerrorthrown");
						}

						// Try here is not necessary since it has been used 20 rows above
						ids_I = tmpNewSpecies; //substrate vector allSpecies position

						// CHECK WHETHER THE NEW PRODUCT IS ALREADY PRESENT AND IF IT ISN'T REACTION IS STORED. IF IT IS, ITS ID BECAME THE VALUE TO STORE INTO THE REACTIONS,
						// OTHERWISE A NEW ID WILL BE CREATED
						ids_II = returnPosSpeciesAlreadyPresent(product_I);
						if(product_I != product_II) // IF PRODUCTS ARE DIFFERENT
						{
							ids_III = returnPosSpeciesAlreadyPresent(product_II);
							// IF BOTH PRODUCTS ARE NEW THE ID OF HTE SECOND ONE IN INCREMENTED BY 1
							if(ids_II == ids_III)
									ids_III++;
						}else{
							ids_III = ids_II;
						}

						// ONE SPECIES CAN NOT CATALYZE BOTH FORWARD AND BACKWARD REACTIONS
						if(toProceed)
						{
							if(!notInverseReactionAlreadyCatalyzed(CLEAVAGE, ids_I, ids_II))
							{
								validReactionFlag = false; // BACKWARD REACTION IS ALREADY CATALIZED!!!
							}else{
							  if(ids_II == (acs_longInt)allSpecies.size()) //IF THE POSITION RETURNED IS EQUAL TO THE SIZE OF vector ALLSPECIES, IT MEANS THAT THE SPECIES IS NOT PRESENT
								{
										if(tmp_RndDoubleGen() < ratioSpeciesEnergizable){tmpEnergizable=ENERGIZABLE;}else{tmpEnergizable=NOTENERGIZABLE;}
										allSpecies.push_back(species(ids_II, product_I, (acs_longInt)0, 0, 0, 0, 0, 0, volume, 0, tmpEnergizable, influx_rate, maxLOut));
								}
							  if(ids_III == (acs_longInt)allSpecies.size()) //IF THE POSITION RETURNED IS EQUAL TO THE SIZE OF vector ALLSPECIES, IT MEANS THAT THE SPECIES IS NOT PRESENT
								{
										if(tmp_RndDoubleGen() < ratioSpeciesEnergizable){tmpEnergizable=ENERGIZABLE;}else{tmpEnergizable=NOTENERGIZABLE;}
										allSpecies.push_back(species(ids_III, product_II, (acs_longInt)0, 0, 0, 0, 0, 0, volume, 0, tmpEnergizable, influx_rate, maxLOut));
								}
								validReactionFlag = true;
							}
						}else{
								validReactionFlag = false;
						}
					}

				//------------------------------ CONDENSATION REACTION -------------------------------------------------|

				}else{
					acs_longInt substratePosID_I;
					acs_longInt substratePosID_II;
					string product;

					if((acs_longInt)tmp_AlreadyEvaluatedSpeciesVector.size() > 0) // IF THERE ARE SPECIES
					{
						// Select if the new species will be the first or the second substrate
						if(tmp_RndDoubleGen() < 0.5)
						{
							substratePosID_I = tmpNewSpecies;
							substratePosID_II = tmp_AlreadyEvaluatedSpeciesVector.at(returnUniformSelection_LONG_IdFromVector(tmp_AlreadyEvaluatedSpeciesVector, tmp_RndDoubleGen));
							tmpCpxTarget = 1;
						}else{
							substratePosID_II = tmpNewSpecies;
							substratePosID_I = tmp_AlreadyEvaluatedSpeciesVector.at(returnUniformSelection_LONG_IdFromVector(tmp_AlreadyEvaluatedSpeciesVector, tmp_RndDoubleGen));
							tmpCpxTarget = 2;
						}

						//PRODUCT CREATION
						try{
							product = allSpecies.at(substratePosID_I).getSequence() + allSpecies.at(substratePosID_II).getSequence();
						}catch(exception&e){
							cout << "allSpecies.at(speciesAvailableForReactions.at(substratePosID_I)).getSequence() + allSpecies.at(speciesAvailableForReactions.at(substratePosID_II)).getSequence();" << endl;
							cout << "Vectorsize "<<allSpecies.size()<<" - position " << substratePosID_I << endl;
							cerr << "exceptioncaught in if(substrateLength > 1):"<<e.what()<<endl;
							ExitWithError("UpdateReactions","exceptionerrorthrown");
						}

						if(debugLevel >= SMALL_DEBUG)
							cout << "\t\t\t\t|- Product: " << product << endl;

						ids_I = returnPosSpeciesAlreadyPresent(product);

						// If the species is one of the initial set
						if(tmpIDtoUpdate <= ((acs_int)firingDisk.size() - 1)){ if(tmpIDtoUpdate == ids_I) toProceed = false;}

						if(toProceed)
						{
							// ONE SPECIES CAN NOT CATALYZE BOTH ONE REACTION AND THE INVERTED ONE
							try{
								ids_II = substratePosID_I;
							}catch(exception&e){
								cout << "speciesAvailableForReactions.at(substratePosID_I)" << endl;
								cout << "Vectorsize "<<tmp_AlreadyEvaluatedSpeciesVector.size()<<" - position "<<substratePosID_I<<endl;
								cerr << "ids_II = speciesAvailableForReactions.at(substratePosID_I);"<<e.what()<<endl;
								ExitWithError("UpdateReactions","exceptionerrorthrown");
							}

							try{
								ids_III = substratePosID_II;
							}catch(exception&e){
								cout << "speciesAvailableForReactions.at(substratePosID_II)" << endl;
								cout << "Vectorsize "<<tmp_AlreadyEvaluatedSpeciesVector.size()<<" - position "<<substratePosID_II<<endl;
								cerr << "ids_II = speciesAvailableForReactions.at(substratePosID_II);"<<e.what()<<endl;
								ExitWithError("UpdateReactions","exceptionerrorthrown");
							}

							if(!notInverseReactionAlreadyCatalyzed(CONDENSATION, ids_I, ids_II))
							{
								validReactionFlag = false;
								if(debugLevel == SMALL_DEBUG)
									 cout << "\t\t\tFound clavage reverse reaction - " << validReactionFlag <<  " -  tentativi - " << rctCreationAttemptCounter << endl;
							}else{
								if(ids_I == (acs_longInt)allSpecies.size()) //IF THE POSITION RETURNED IS EQUAL TO THE SIZE OF VECOTR ALLSPECIES, IT MEANS THAT THE SPECIES IS NOT PRESENT
								{
										if(tmp_RndDoubleGen() < ratioSpeciesEnergizable){tmpEnergizable=ENERGIZABLE;}else{tmpEnergizable=NOTENERGIZABLE;}
										allSpecies.push_back(species(ids_I, product, (acs_longInt)0, 0, 0, 0, 0, 0, volume, 0, tmpEnergizable, influx_rate, maxLOut));
								}
								validReactionFlag = true;
							}
						}else{ // if(toProceed)
							validReactionFlag = false;
						}
					}//if((acs_longInt)tmp_AlreadyEvaluatedSpeciesVector.size() > 0) // IF THERE ARE SPECIES
				}

				if(debugLevel == SMALL_DEBUG)
					cout << "\t\t\tReaction validity: " << validReactionFlag << " - attempts " << rctCreationAttemptCounter << endl;

				// IF ALL IT'S OK REACTION AND CATALYSIS ARE STORED
				if(validReactionFlag == true)
				{
					// SINCE THE SPECIES HAS BEEN CREATED, NOW I HAVE TO CHECK WHETHER THE REACTION IS ALREADY PRESENT.
					if(getNumberOfReactions() > 0)
					{
						id_reaction = returnPosReactionAlreadyPresent(tmpRctType, ids_I, ids_II, ids_III);
					}else{
						id_reaction = 0;
					}

					if(id_reaction == getNumberOfReactions())
					{
						// IF THE ENERGY IS CONSIDERED
						// According to the energy parameter and the reaction type the energetic type is assigned
						if(nrgBoolFlag == ENERGYBASED)
						{
							if(nrgBooleanFunctions.size() > 0)
							{
								nrgRctBoolFncID = returnSelectionIdFromAWeightProbVectorAlreadyNormalized(nrgBoolFncsProb2BeSelected,tmp_RndDoubleGen);
								try{
									tempEnergyRctType = nrgBooleanFunctions.at(nrgRctBoolFncID);
								}catch(exception&e){
									cout << "tempEnergyRctType = nrgBooleanFunctions.at(nrgRctBoolFncID);" << endl;
									cout << "Vectorsize "<<nrgBooleanFunctions.size()<<" - position "<<nrgRctBoolFncID<<endl;
									cerr << "tempEnergyRctType = nrgBooleanFunctions.at(nrgRctBoolFncID);"<<e.what()<<endl;
									ExitWithError("createReactionsForThisSpecies","exceptionerrorthrown");
								}
							}else{
								// If the vector contains only 1 element, the first one is selected
								tempEnergyRctType = nrgBooleanFunctions.front();
							}
						}
						// IF IT IS NOT PRESENT THE VECTOR ALLREACTIONS AND THE VECTOR ALLCATALYSIS ARE UPDATED
                        if(spontRct == NOSPONTANEOUS){tmpK_spontanoues = NOSPONTANEOUS;}
                        else{
                        	if((tmpRctType == CONDENSATION) || (tmpRctType == ENDO_CONDENSATION)){tmpK_spontanoues = K_spont_ass;}
                        	else{tmpK_spontanoues = K_spont_diss;}
                        }
						allReactions.push_back(reactions(id_reaction, tmpRctType, ids_I, ids_II, ids_III, 0, tempEnergyRctType, tmpK_spontanoues));
						if(reverseReactions == true)
						{
							if(tmpRctType == CONDENSATION)
							{
								tmpK_ass = K_ass;
								tmpK_diss = K_diss / revRctRatio;
								tmpK_cpx = K_cpx;
							}else{ //Cleavage reaction
								tmpK_ass = K_ass / revRctRatio;
								tmpK_diss = K_diss;
								tmpK_cpx = K_cpx / revRctRatio;
							}
						}else{ // reverseReactions == FALSE
							tmpK_ass = K_ass;
							tmpK_diss = K_diss;
							tmpK_cpx = K_cpx;
						}
						// IF camplex target must be defined as 1
						if(allSpecies.at(tmpIDtoUpdate).getSequenceLength() <= nonCatalyticMaxLength) ExitWithError("createReactionsForThisSpecies","ERROR: The catalyst is smaller than the minumum");
						if((tmpCpxTarget < 1) || (tmpCpxTarget > 2))
						{
							cout << "id reaction: " << id_reaction << " - Reaction type: " << allReactions.at(id_reaction).getType() << endl;
							ExitWithError("createReactionsForThisSpecies:r1154","WRONG COMPLEX TARGET");
						}
						allCatalysis.push_back(catalysis(getNumberOfCatalysis(), tmpIDtoUpdate, id_reaction, 0, tmpK_ass, tmpK_diss, tmpK_cpx, tmpCpxTarget));

					}else{ // OTHERWISE ONLY THE VECTOR ALLCATALYSIS IS UPDATED
						if(reverseReactions == true)
						{
							if(tmpRctType == CONDENSATION)
							{
								tmpK_ass = K_ass;
								tmpK_diss = K_diss / revRctRatio;
								tmpK_cpx = K_cpx;
							}else{ //Cleavage reaction
								tmpK_ass = K_ass / revRctRatio;
								tmpK_diss = K_diss;
								tmpK_cpx = K_cpx / revRctRatio;
							}
						}else{ // reverseReactions == FALSE
							tmpK_ass = K_ass;
							tmpK_diss = K_diss;
							tmpK_cpx = K_cpx;
						}
						//Check if this reaction is not already catalysed by this species
						if(!checkIfTheReactionIsAlreadyCatalyzedByThisSpecies(tmpIDtoUpdate, id_reaction))
						{
							if(allSpecies.at(tmpIDtoUpdate).getSequenceLength() <= nonCatalyticMaxLength) ExitWithError("createReactionsForThisSpecies","ERROR: The catalyst is smaller than the minumum");
							if((tmpCpxTarget < 1) || (tmpCpxTarget > 2)) ExitWithError("createReactionsForThisSpecies:r1179","WRONG COMPLEX TARGET");
							allCatalysis.push_back(catalysis(getNumberOfCatalysis(), tmpIDtoUpdate, id_reaction, 0, tmpK_ass, tmpK_diss, tmpK_cpx, tmpCpxTarget));
						}
					}
					if(debugLevel >= SMALL_DEBUG)
					{
						cout << "\t\t\t|- old species ID: " << allSpecies.at(tmpIDtoUpdate).getID() << " - Sequence: " << allSpecies.at(tmpIDtoUpdate).getSequence() << endl;
						cout << "\t\t\t\t have been updated with new catalysis: " << getNumberOfCatalysis() << " affecting reaction " << id_reaction << "(" << tmpRctType << ")"
								<< " involving species " << ids_I << " (" << allSpecies.at(ids_I).getSequence() << ") - "
								<< ids_II  << " (" << allSpecies.at(ids_II).getSequence() << ") - "
								<< ids_III << " (" << allSpecies.at(ids_III).getSequence() << ") - complex target: " << tmpCpxTarget << endl;
					}
				}else{ // if(validReactionFlag == true)
						rctCreationAttemptCounter++;
				}// if(validReactionFlag)
		} //while(!validReactionFlag)

	if(debugLevel == FINDERRORDURINGRUNTIME) cout << "\t\t\t\tenvironment::createReactionsForThisSpecies end" << endl;
	return true;
}

/**
 Select reaction type (cleavage or condensation) according to the cleavage probability
 @version 1.0
 @param acs_int tmpTotalNumberOfReactions Total number of conceivable reactions
 @param acs_double tmpRctsProb reaction probability
 */
acs_int environment::selectWhetherCleavageOrCond(MTRand& tmp__RndDoubleGen)
{
	if(debugLevel == FINDERRORDURINGRUNTIME) cout << "\t\t\t\t\tenvironment::selectWhetherCleavageOrCond start" << endl;
	
	acs_int tmpReactionType = CONDENSATION; //set condensation as default value
	acs_double rndCleavProm = tmp__RndDoubleGen();	
	if(rndCleavProm < cleavageProbability)
	{
		tmpReactionType = CLEAVAGE;
	}	
	
	if(debugLevel == FINDERRORDURINGRUNTIME) cout << "\t\t\t\t\tenvironment::selectWhetherCleavageOrCond end" << endl;

	return tmpReactionType;
}

/**
 Initial molecule population creation from file using standard C++ libraries
 @version 1.0
 @param string tmpSpeciesFilePath file path
 @date 20130702
 */
bool environment::createInitialMoleculesPopulationFromFileSTD(string tmpSpeciesFilePath)
{
    if(debugLevel == FINDERRORDURINGRUNTIME) cout << "environment::createInitialMoleculesPopulationFromFileSTD start" << endl;
    // SPECIES FILE PATH CREATION
    bool maxNonCatIdDef = false;
    string SpeciesFilePath = tmpSpeciesFilePath + "_acsspecies.csv";
    cout << "|- Loading Species from " << SpeciesFilePath;
    ifstream myfile;
    try{
    	myfile.open(SpeciesFilePath.c_str());
    }catch(exception&e){
		cerr << "|- ERROR!!! PROBLEM WITH FILE " << SpeciesFilePath << e.what() << endl;
		ExitWithError("createInitialMoleculesPopulationFromFile","exceptionerrorthrown");
    }
    string strID, strCod, strConc, strDiff, strPrec, strK_cpx, strCpxBin, strEval, strAge, strReb, strCatID, strCpxID, strPho, strChar, strLock;
    while (myfile.good())
    {
        getline(myfile, strID, '\t');
        getline(myfile, strCod, '\t');
        getline(myfile, strConc, '\t');
        getline(myfile, strDiff, '\t');
        getline(myfile, strPrec, '\t');
        getline(myfile, strK_cpx, '\t');
        getline(myfile, strCpxBin, '\t');
        getline(myfile, strEval, '\t');
        getline(myfile, strAge, '\t');
        getline(myfile, strReb, '\t');
        getline(myfile, strCatID, '\t');
        getline(myfile, strCpxID, '\t');
        getline(myfile, strPho, '\t');
        getline(myfile, strChar, '\t');
        getline(myfile, strLock, '\n');

        if(strID.find("\n") != 0 && (strID.size() > 0) && (strID.find(" ") != 0))
        {
			allSpecies.push_back(species((acs_longInt)atol(strID.c_str()), strCod, (acs_double)atof(strConc.c_str()),
										 (acs_double)atof(strDiff.c_str()),(acs_int)atoi(strPrec.c_str()),
										 (acs_double)atof(strK_cpx.c_str()), (acs_int)atoi(strCpxBin.c_str()),
										 (acs_int)atoi(strEval.c_str()), (acs_double)atof(strAge.c_str()),atoi(strReb.c_str()), volume,
										 (acs_longInt)atol(strCatID.c_str()), (acs_longInt)atol(strCpxID.c_str()),
										 (acs_double)atof(strPho.c_str()), (acs_double)atof(strChar.c_str()),
										 (acs_int)atoi(strLock.c_str()), influx_rate, maxLOut));
			// Define the max non catalytic max ID
			if(!maxNonCatIdDef){if(strCod.size() > nonCatalyticMaxLength){nonCatalyticLastID = (acs_longInt)atol(strID.c_str()) - 1; maxNonCatIdDef=true;}}
			try{
				if(allSpecies.at((acs_longInt)atoi(strID.c_str())).getComplexCutPnt() == 0)
				{
					if(allSpecies.at((acs_longInt)atoi(strID.c_str())).getAmount() > 0)
						numberOfSpecies++;

					numberOfMolecules += allSpecies.at((acs_longInt)atoi(strID.c_str())).getAmount();
					if(((acs_longInt)atoi(strID.c_str()) > lastFiringDiskSpeciesID) &&
					   (allSpecies.at((acs_longInt)atoi(strID.c_str())).getAmount() > 0))
					{
						numberOfNewSpecies++;
						numberOfNewMolecules += allSpecies.at((acs_longInt)atoi(strID.c_str())).getAmount();
					}
				}else{ // If the species is a complex
					incNumberOfCpx();
					numberOfCpxMols += allSpecies.at((acs_longInt)atoi(strID.c_str())).getAmount();
				}
			}catch(exception&e){
				cout << "allSpecies.at((acs_longInt)strLineSpletted[0].toInt()).getComplexCutPnt()..." << endl;
				cout << "Vectorsize " << allSpecies.size() << " - position "<< (acs_longInt)atoi(strID.c_str()) << endl;
				cerr << "allSpecies.at((acs_longInt)strLineSpletted[0].toInt()).getComplexCutPnt().."<<e.what()<<endl;
				ExitWithError("createInitialMoleculesPopulationFromFile","exceptionerrorthrown");
			}
        }
    }
    myfile.close();
    // FIRING DISK STORING, NUTRIENTS CREATION AND PROBABILITY VECTOR FOR A NUTRIENT TO BE SELECTED IS CREATED
    for(acs_int singleSpecies = 0; singleSpecies < (acs_int)allSpecies.size(); singleSpecies++)
    {
        if(allSpecies.at(singleSpecies).getID() <= getLastFiringDiskSpeciesID())
            firingDisk.push_back(allSpecies.at(singleSpecies));
        if(debugLevel == SMALL_DEBUG)
        	cout << allSpecies.at(singleSpecies).getID() << endl;
    }
    cout << " Done" << endl;
    if(debugLevel == FINDERRORDURINGRUNTIME) cout << "environment::createInitialMoleculesPopulationFromFile end" << endl;

    return true;
}//eof createInitialPopulationFromFileSTD



/**
 Initial molecule population creation from file
 @version 1.0
 @param QString tmpSpeciesFilePath file path
 */
/*bool environment::createInitialMoleculesPopulationFromFile(QString tmpSpeciesFilePath)
{
	if(debugLevel == FINDERRORDURINGRUNTIME) cout << "environment::createInitialMoleculesPopulationFromFile start" << endl;

	// CONFIGURATION FILE PATH CREATION
	QString SpeciesFilePath = tmpSpeciesFilePath + "_acsspecies.csv";
	QTextStream prompt(stdout);
	prompt << "Loading species sequences from " << SpeciesFilePath << endl;
	//OPEN FILE
	QFile fid(SpeciesFilePath);
	if(!fid.open(QIODevice::ReadOnly | QIODevice::Text))
		ExitWithError("createInitialMoleculesPopulationFromFile", "File has not been found");
	
	QTextStream in(&fid);
	while(!in.atEnd())
	{
            QString strLine = in.readLine();
            if((strLine.contains("#") == false) && (strLine.isEmpty() == false) && (strLine.isNull() == false))
            {
                QStringList strLineSpletted = strLine.split("\t");
                allSpecies.push_back(species((acs_longInt)strLineSpletted[0].toInt(), strLineSpletted[1].QString::toStdString(),
                                             (acs_double)strLineSpletted[2].toDouble(), (acs_double)strLineSpletted[3].toDouble(),
                                             (acs_double)strLineSpletted[4].toDouble(), (acs_double)strLineSpletted[5].toDouble(),
                                             (acs_longInt)strLineSpletted[6].toInt(), strLineSpletted[7].toInt(),
                                             (acs_double)strLineSpletted[8].toDouble(), strLineSpletted[9].toInt(), volume,
                                             (acs_double)strLineSpletted[10].toDouble(), (acs_double)strLineSpletted[11].toDouble(),
                                             (acs_double)strLineSpletted[12].toDouble(), (acs_double)strLineSpletted[13].toDouble(),
                                             (acs_int)strLineSpletted[14].toInt(), influx_rate, maxLOut));
                try{
                    if(allSpecies.at((acs_longInt)strLineSpletted[0].toInt()).getComplexCutPnt() == 0)
                    {
                        if(allSpecies.at((acs_longInt)strLineSpletted[0].toInt()).getAmount() > 0)
                            numberOfSpecies++;

                        numberOfMolecules += allSpecies.at((acs_longInt)strLineSpletted[0].toInt()).getAmount();
                        if(((acs_longInt)strLineSpletted[0].toInt() > lastFiringDiskSpeciesID) &&
                           (allSpecies.at((acs_longInt)strLineSpletted[0].toInt()).getAmount() > 0))
                        {
                            numberOfNewSpecies++;
                            numberOfNewMolecules += allSpecies.at((acs_longInt)strLineSpletted[0].toInt()).getAmount();
                        }
                    }else{ // If the species is a complex
                        incNumberOfCpx();
                        numberOfCpxMols += allSpecies.at((acs_longInt)strLineSpletted[0].toInt()).getAmount();
                    }
                }catch(exception&e){
                    cout << "allSpecies.at((acs_longInt)strLineSpletted[0].toInt()).getComplexCutPnt()..." << endl;
                    cout << "Vectorsize " << allSpecies.size() << " - position "<< (acs_longInt)strLineSpletted[0].toInt() << endl;
                    cerr << "allSpecies.at((acs_longInt)strLineSpletted[0].toInt()).getComplexCutPnt().."<<e.what()<<endl;
                    ExitWithError("createInitialMoleculesPopulationFromFile","exceptionerrorthrown");
                }
            }
        }
	
	fid.close();
	// FIRING DISK STORING, NUTRIENTS CREATION AND PROBABILITY VECTOR FOR A NUTRIENT TO BE SELECTED IS CREATED
	for(acs_int singleSpecies = 0; singleSpecies < (acs_int)allSpecies.size(); singleSpecies++)
	{
		if(allSpecies.at(singleSpecies).getID() <= getLastFiringDiskSpeciesID())
			firingDisk.push_back(allSpecies.at(singleSpecies));
	}	
	if(debugLevel == FINDERRORDURINGRUNTIME) cout << "environment::createInitialMoleculesPopulationFromFile end" << endl;
	
	return true;
}//eof createInitialPopulationFromFile*/

/**
 Initial molecule population creation. Species are uploaed from a SPECIFIC file created using actual generation and simuation
 @version 1.0
 @param QString tmpSpeciesFilePath file path
 */
/*bool environment::createInitialMoleculesPopulationFromSpecificFile(QString tmpSpeciesFilePath, acs_int tmpActGEN, acs_int tmpActSIM)
{
	if(debugLevel == FINDERRORDURINGRUNTIME) cout << "environment::createInitialMoleculesPopulationFromSpecificFile start" << endl;
	
	QString strCurrentGen;
	QString strCurrentSim;
	QString strCurrentStep;
	QString strZeroGenBefore = zeroBeforeStringNumber(nGEN, tmpActGEN);
	QString strZeroSimBefore = zeroBeforeStringNumber(pow(double(nSIM), double(tmpActGEN)), tmpActSIM);
	QString strZeroStepBefore = zeroBeforeStringNumber(nReactions, 0);
	// CONFIGURATION FILE PATH CREATION
	QString SpeciesFilePath = tmpSpeciesFilePath + "species_" + strZeroGenBefore + strCurrentGen.setNum(tmpActGEN) + 
															"_" +  strZeroSimBefore + strCurrentSim.setNum(tmpActSIM) + 
															"_" +  strZeroStepBefore + strCurrentStep.setNum(0)	+ ".csv";
	QTextStream prompt(stdout);
	if(debugLevel >= MINIMAL_PROMPT)
	{
		prompt << endl << "\t\t|- Loading species structure from " << SpeciesFilePath << endl;
	}
	//OPEN FILE
	QFile fid(SpeciesFilePath);
	if(!fid.open(QIODevice::ReadOnly | QIODevice::Text))
		ExitWithError("createInitialMoleculesPopulationFromSpecificFile", "File has not been found");
	
	QTextStream in(&fid);
	while(!in.atEnd())
	{
		QString strLine = in.readLine();
		if((strLine.contains("#") == false) && (strLine.isEmpty() == false) && (strLine.isNull() == false))
		{
			QStringList strLineSpletted = strLine.split("\t");
			allSpecies.push_back(species((acs_longInt)strLineSpletted[0].toInt(), strLineSpletted[1].QString::toStdString(), 
                                                     (acs_double)strLineSpletted[2].toDouble(), (acs_double)strLineSpletted[3].toDouble(),
                                                     (acs_double)strLineSpletted[4].toDouble(), (acs_double)strLineSpletted[5].toDouble(),
                                                     (acs_longInt)strLineSpletted[6].toInt(), strLineSpletted[7].toInt(),
                                                     (acs_double)strLineSpletted[8].toDouble(), strLineSpletted[9].toInt(), volume,
                                                     (acs_double)strLineSpletted[10].toDouble(), (acs_double)strLineSpletted[11].toDouble(),
                                                     (acs_double)strLineSpletted[12].toDouble(), (acs_double)strLineSpletted[13].toDouble(),
                                                     (acs_int)strLineSpletted[14].toInt(), influx_rate, maxLOut));
			if(allSpecies.at((acs_longInt)strLineSpletted[0].toInt()).getComplexCutPnt() == 0)
			{
                if(allSpecies.at((acs_longInt)strLineSpletted[0].toInt()).getAmount() > 0) numberOfSpecies++;
				numberOfMolecules += allSpecies.at((acs_longInt)strLineSpletted[0].toInt()).getAmount();
                if(((acs_longInt)strLineSpletted[0].toInt() > lastFiringDiskSpeciesID) && (allSpecies.at((acs_longInt)strLineSpletted[0].toInt()).getAmount() > 0))
				{
					numberOfNewSpecies++;
					numberOfNewMolecules += allSpecies.at((acs_longInt)strLineSpletted[0].toInt()).getAmount();
				}
			}else{ // If the species is a complex
				incNumberOfCpx();
				numberOfCpxMols += allSpecies.at((acs_longInt)strLineSpletted[0].toInt()).getAmount();
			}
		}
	}	
	
	fid.close();
	
	if(debugLevel == FINDERRORDURINGRUNTIME) cout << "environment::createInitialMoleculesPopulationFromSpecificFile end" << endl;
	
	return true;
}//eof createInitialPopulationFromFile*/

/**
 Initial molecule population creation. Species are uploaed from a SPECIFIC file created using actual generation and simuation (standard C++ libraries)
 @version 1.0
 @param QString tmpSpeciesFilePath file path
 */
bool environment::createInitialMoleculesPopulationFromSpecificFileSTD(string tmpSpeciesFilePath, acs_int tmpActGEN, acs_int tmpActSIM)
{
    if(debugLevel == FINDERRORDURINGRUNTIME) cout << "environment::createInitialMoleculesPopulationFromSpecificFileSTD start" << endl;

    stringstream strCurrentGen;
    stringstream strCurrentSim;
    stringstream strCurrentStep;
    bool maxNonCatIdDef = false;
    string strZeroGenBefore = zeroBeforeStringNumberSTD(nGEN, tmpActGEN);
    string strZeroSimBefore = zeroBeforeStringNumberSTD(pow(double(nSIM), double(tmpActGEN)), tmpActSIM);
    string strZeroStepBefore = zeroBeforeStringNumberSTD(nReactions, 0);

    strCurrentGen << tmpActGEN;
    strCurrentSim << tmpActSIM;
    strCurrentStep << 0;
    string SpeciesFilePath = tmpSpeciesFilePath + "species_" + strZeroGenBefore + strCurrentGen.str() +
            "_" +  strZeroSimBefore + strCurrentSim.str() +
            "_" +  strZeroStepBefore + strCurrentStep.str()	+ ".csv";

    ifstream myfile;
    myfile.open(SpeciesFilePath.c_str());
    string strID, strCod, strConc, strDiff, strPrec, strK_cpx, strCpxBin, strEval, strAge, strReb, strCatID, strCpxID, strPho, strChar, strLock;
    while (myfile.good())
    {
        getline(myfile, strID, '\t');
        getline(myfile, strCod, '\t');
        getline(myfile, strConc, '\t');
        getline(myfile, strDiff, '\t');
        getline(myfile, strPrec, '\t');
        getline(myfile, strK_cpx, '\t');
        getline(myfile, strCpxBin, '\t');
        getline(myfile, strEval, '\t');
        getline(myfile, strAge, '\t');
        getline(myfile, strReb, '\t');
        getline(myfile, strCatID, '\t');
        getline(myfile, strCpxID, '\t');
        getline(myfile, strPho, '\t');
        getline(myfile, strChar, '\t');
        getline(myfile, strLock, '\n');

        if(strID.find("\n") != 0 && (strID.size() > 0) && (strID.find(" ") != 0))
        {
			allSpecies.push_back(species((acs_longInt)atol(strID.c_str()), strCod, (acs_double)atof(strConc.c_str()),
										 (acs_double)atof(strDiff.c_str()),(acs_int)atoi(strPrec.c_str()),
										 (acs_double)atof(strK_cpx.c_str()), (acs_int)atoi(strCpxBin.c_str()),
										 (acs_int)atoi(strEval.c_str()), (acs_double)atof(strAge.c_str()),atoi(strReb.c_str()), volume,
										 (acs_longInt)atol(strCatID.c_str()), (acs_longInt)atol(strCpxID.c_str()),
										 (acs_double)atof(strPho.c_str()), (acs_double)atof(strChar.c_str()),
										 (acs_int)atoi(strLock.c_str()), influx_rate, maxLOut));
			if(!maxNonCatIdDef){if(strCod.size() > nonCatalyticMaxLength){nonCatalyticLastID = (acs_longInt)atol(strID.c_str()) - 1; maxNonCatIdDef=true;}}
			try{
				if(allSpecies.at((acs_longInt)atoi(strID.c_str())).getComplexCutPnt() == 0)
				{
					if(allSpecies.at((acs_longInt)atoi(strID.c_str())).getAmount() > 0)
						numberOfSpecies++;

					numberOfMolecules += allSpecies.at((acs_longInt)atoi(strID.c_str())).getAmount();
					if(((acs_longInt)atoi(strID.c_str()) > lastFiringDiskSpeciesID) &&
					   (allSpecies.at((acs_longInt)atoi(strID.c_str())).getAmount() > 0))
					{
						numberOfNewSpecies++;
						numberOfNewMolecules += allSpecies.at((acs_longInt)atoi(strID.c_str())).getAmount();
					}
				}else{ // If the species is a complex
					incNumberOfCpx();
					numberOfCpxMols += allSpecies.at((acs_longInt)atoi(strID.c_str())).getAmount();
				}
			}catch(exception&e){
				cout << "allSpecies.at((acs_longInt)strLineSpletted[0].toInt()).getComplexCutPnt()..." << endl;
				cout << "Vectorsize " << allSpecies.size() << " - position "<< (acs_longInt)atoi(strID.c_str()) << endl;
				cerr << "allSpecies.at((acs_longInt)strLineSpletted[0].toInt()).getComplexCutPnt().."<<e.what()<<endl;
				ExitWithError("createInitialMoleculesPopulationFromFile","exceptionerrorthrown");
			}
        }
    }
    myfile.close();
    // FIRING DISK STORING, NUTRIENTS CREATION AND PROBABILITY VECTOR FOR A NUTRIENT TO BE SELECTED IS CREATED
    for(acs_int singleSpecies = 0; singleSpecies < (acs_int)allSpecies.size(); singleSpecies++)
    {
        if(allSpecies.at(singleSpecies).getID() <= getLastFiringDiskSpeciesID())
            firingDisk.push_back(allSpecies.at(singleSpecies));
    }

    if(debugLevel == FINDERRORDURINGRUNTIME) cout << "environment::createInitialMoleculesPopulationFromSpecificFileSTD end" << endl;

    return true;
}//eof createInitialMoleculesPopulationFromSpecificFileSTD

/**
 Create influx layer from file C++ libraries
 @version 1.0
 @param string tmpInfluxFilePath file path
 @date 20130702
 */
bool environment::createInfluxLayersFromFileSTD(string tmpInfluxFilePath)
{
    if(debugLevel == FINDERRORDURINGRUNTIME) cout << "environment::createInfluxLayersFromFileSTD start" << endl;
    // SPECIES FILE PATH CREATION
    string FilePath = tmpInfluxFilePath + "_acsinflux.csv";
    ifstream myfile;
    myfile.open(FilePath.c_str());
    string strFeedID, strFeedProb;
    while (myfile.good())
    {
        getline(myfile, strFeedID, '\t');
        getline(myfile, strFeedProb, '\n');

        if(strFeedID.find("\n") != 0 && (strFeedID.size() > 0) && (strFeedID.find(" ") != 0))
        {
        	nutrientsForInflux.push_back((acs_int)atoi(strFeedID.c_str()));
        	nutrientsProb2BeSelected.push_back((acs_double)atof(strFeedProb.c_str()));
        }

    }

    myfile.close();
    if(debugLevel == FINDERRORDURINGRUNTIME) cout << "environment::createInfluxLayersFromFileSTD end" << endl;
    return true;
}//eof createInfluxLayersFromFileSTD

/**
 load energy boolean function (in decimal format) - Standard C++
 @version 1.0
 @param string tmpBoolNrgFilePath file path
 @date 20130702
 */
bool environment::createNrgBooleanFunctionsFromFileSTD(string tmpBoolNrgFilePath)
{
    if(debugLevel == FINDERRORDURINGRUNTIME) cout << "environment::createNrgBooleanFunctionsFromFileSTD start" << endl;
    // SPECIES FILE PATH CREATION
    string FilePath = tmpBoolNrgFilePath + "_acsnrgbooleanfunctions.csv";
    ifstream myfile;
    myfile.open(FilePath.c_str());
    string strBoolNrgID, strBoolNrgProb;
    while (myfile.good())
    {
        getline(myfile, strBoolNrgID, '\t');
        getline(myfile, strBoolNrgProb, '\n');

        if(strBoolNrgID.find("\n") != 0 && (strBoolNrgID.size() > 0) && (strBoolNrgID.find(" ") != 0))
        {
        	nrgBooleanFunctions.push_back((acs_int)atoi(strBoolNrgID.c_str()));
        	nrgBoolFncsProb2BeSelected.push_back((acs_double)atof(strBoolNrgProb.c_str()));
        }
    }

    myfile.close();
    if(debugLevel == FINDERRORDURINGRUNTIME) cout << "environment::createNrgBooleanFunctionsFromFileSTD end" << endl;
    return true;
}//eof createNrgBooleanFunctionsFromFileSTD

/**
 Reactions from file using standard C++ libraries
 @version 1.0
 @param string tmpSpeciesFilePath file path
 @date 20130702
 */
bool environment::createInitialReactionsLayerFromFileSTD(string tmpSpeciesFilePath)
{
    if(debugLevel == FINDERRORDURINGRUNTIME) cout << "environment::createInitialReactionsLayerFromFileSTD start" << endl;
    // SPECIES FILE PATH CREATION
    cout << "|- Loading Reactions...";
    string FilePath = tmpSpeciesFilePath + "_acsreactions.csv";
    ifstream myfile;
    myfile.open(FilePath.c_str());
    string strID, strType, strMolsI, strMolsII, strMolsIII, strCnt, strNrg, strKspont;
    while (myfile.good())
    {

        getline(myfile, strID, '\t');
        getline(myfile, strType, '\t');
        getline(myfile, strMolsI, '\t');
        getline(myfile, strMolsII, '\t');
        getline(myfile, strMolsIII, '\t');
        getline(myfile, strCnt, '\t');
        getline(myfile, strNrg, '\t');
        getline(myfile, strKspont, '\n');
        if((strID.find("\n") != 0) && (strID.size() > 0) && (strID.find(" ") != 0))
            allReactions.push_back(reactions((acs_longInt)atol(strID.c_str()), (acs_int)atoi(strType.c_str()), (acs_longInt)atol(strMolsI.c_str()),
                                         (acs_longInt)atol(strMolsII.c_str()), (acs_longInt)atol(strMolsIII.c_str()), (acs_longInt)atol(strCnt.c_str()),
                                         (acs_double)atof(strNrg.c_str()), (acs_double)atof(strKspont.c_str())));

    }

    myfile.close();

    if(debugLevel == FINDERRORDURINGRUNTIME) cout << "environment::createInitialReactionsLayerFromFileSTD end" << endl;
    return true;
}//eof createInitialReactionsLayerFromFileSTD

/**
 Reactions from file using standard C++ libraries
 @version 1.0
 @param string tmpSpeciesFilePath file path
 @date 20130702
 */
bool environment::createInitialReactionsLayerFromSpecificFileSTD(string tmpReactionsFilePath, acs_int tmpActGEN, acs_int tmpActSIM)
{
    if(debugLevel == FINDERRORDURINGRUNTIME) cout << "environment::createInitialReactionsLayerFromSpecificFileSTD start" << endl;
    // CONFIGURATION FILE PATH CREATION
    stringstream strCurrentGen;
    stringstream strCurrentSim;
    stringstream strCurrentStep;
    string strZeroGenBefore = zeroBeforeStringNumberSTD(nGEN, tmpActGEN);
    string strZeroSimBefore = zeroBeforeStringNumberSTD(pow(double(nSIM), double(tmpActGEN)), tmpActSIM);
    string strZeroStepBefore = zeroBeforeStringNumberSTD(nReactions, 0);

    strCurrentGen << tmpActGEN;
    strCurrentSim << tmpActSIM;
    strCurrentStep << 0;
    string ReactionsFilePath = tmpReactionsFilePath + "reactions_" + strZeroGenBefore + strCurrentGen.str() +
            "_" +  strZeroSimBefore + strCurrentSim.str() +
            "_" +  strZeroStepBefore + strCurrentStep.str()	+ ".csv";

    ifstream myfile;
    myfile.open(ReactionsFilePath.c_str());
    string strID, strType, strMolsI, strMolsII, strMolsIII, strCnt, strNrg, strKspont;
    while (myfile.good())
    {
        getline(myfile, strID, '\t');
        getline(myfile, strType, '\t');
        getline(myfile, strMolsI, '\t');
        getline(myfile, strMolsII, '\t');
        getline(myfile, strMolsIII, '\t');
        getline(myfile, strCnt, '\t');
        getline(myfile, strNrg, '\t');
        getline(myfile, strKspont, '\n');

        if(strID.find("\n") != 0 && (strID.size() > 0) && (strID.find(" ") != 0))
        	allReactions.push_back(reactions((acs_longInt)atol(strID.c_str()), (acs_int)atoi(strType.c_str()), (acs_longInt)atol(strMolsI.c_str()),
                                         (acs_longInt)atol(strMolsII.c_str()), (acs_longInt)atol(strMolsIII.c_str()), (acs_longInt)atol(strCnt.c_str()),
                                         (acs_double)atof(strNrg.c_str()), (acs_double)atof(strKspont.c_str())));

    }

    myfile.close();
    if(debugLevel == FINDERRORDURINGRUNTIME) cout << "environment::createInitialReactionsLayerFromSpecificFileSTD end" << endl;
    return true;
}//eof createInitialReactionsLayerFromSpecificFileSTD

/**
 Catalysis from file using standard C++ libraries
 @version 1.0
 @param string tmpSpeciesFilePath file path
 @date 20130702
 */
bool environment::createInitialCatalysisLayerFromFileSTD(string tmpCatalysisFilePath)
{
    if(debugLevel == FINDERRORDURINGRUNTIME) cout << "environment::createInitialCatalysisLayerFromFileSTD start" << endl;
    // SPECIES FILE PATH CREATION
    cout << "|- Loading Catalysis..." << endl;
    string FilePath = tmpCatalysisFilePath + "_acscatalysis.csv";
    ifstream myfile;
    myfile.open(FilePath.c_str());
    string strID, strCatalyst, strReaction, strCounter, strKcond, strKcleav, strKcpx, strCpxTarget;
    while (myfile.good())
    {
        getline(myfile, strID, '\t');
        getline(myfile, strCatalyst, '\t');
        getline(myfile, strReaction, '\t');
        getline(myfile, strCounter, '\t');
        getline(myfile, strKcond, '\t');
        getline(myfile, strKcleav, '\t');
        getline(myfile, strKcpx, '\t');
        getline(myfile, strCpxTarget, '\n');

        if(strID.find("\n") != 0 && (strID.size() > 0) && (strID.find(" ") != 0))
        {
        	allCatalysis.push_back(catalysis((acs_longInt)atol(strID.c_str()), (acs_longInt)atol(strCatalyst.c_str()), (acs_longInt)atol(strReaction.c_str()),
                                         (acs_longInt)atol(strCounter.c_str()), (acs_double)atof(strKcond.c_str()), (acs_double)atof(strKcleav.c_str()),
                                         (acs_double)atof(strKcpx.c_str()), (acs_int)atoi(strCpxTarget.c_str())));
        }

    }

    myfile.close();
    inserSubListInSpecies();
    cout << " Done" << endl;
    if(debugLevel == FINDERRORDURINGRUNTIME) cout << "environment::createInitialCatalysisLayerFromFileSTD end" << endl;
    return true;
}//eof createInitialCatalysisLayerFromFileSTD

/**
 catalysis from file using standard C++ libraries
 @version 1.0
 @param string tmpSpeciesFilePath file path
 @date 20130702
 */
bool environment::createInitialCatalysisLayerFromSpecificFileSTD(string tmpCatalysisFilePath, acs_int tmpActGEN, acs_int tmpActSIM)
{
    if(debugLevel == FINDERRORDURINGRUNTIME) cout << "environment::createInitialCatalysisLayerFromSpecificFileSTD start" << endl;
    // CONFIGURATION FILE PATH CREATION
    stringstream strCurrentGen;
    stringstream strCurrentSim;
    stringstream strCurrentStep;
    string strZeroGenBefore = zeroBeforeStringNumberSTD(nGEN, tmpActGEN);
    string strZeroSimBefore = zeroBeforeStringNumberSTD(pow(double(nSIM), double(tmpActGEN)), tmpActSIM);
    string strZeroStepBefore = zeroBeforeStringNumberSTD(nReactions, 0);

    strCurrentGen << tmpActGEN;
    strCurrentSim << tmpActSIM;
    strCurrentStep << 0;
    string catalysisFilePath = tmpCatalysisFilePath + "catalysis_" + strZeroGenBefore + strCurrentGen.str() +
            "_" +  strZeroSimBefore + strCurrentSim.str() +
            "_" +  strZeroStepBefore + strCurrentStep.str()	+ ".csv";

    ifstream myfile;
    myfile.open(catalysisFilePath.c_str());
    string strID, strCatalyst, strReaction, strCounter, strKcond, strKcleav, strKcpx, strCpxTarget;
    while (myfile.good())
    {
        getline(myfile, strID, '\t');
        getline(myfile, strCatalyst, '\t');
        getline(myfile, strReaction, '\t');
        getline(myfile, strCounter, '\t');
        getline(myfile, strKcond, '\t');
        getline(myfile, strKcleav, '\t');
        getline(myfile, strKcpx, '\t');
        getline(myfile, strCpxTarget, '\n');

        if(strID.find("\n") != 0 && (strID.size() > 0) && (strID.find(" ") != 0))
        	allCatalysis.push_back(catalysis((acs_longInt)atol(strID.c_str()), (acs_longInt)atol(strCatalyst.c_str()), (acs_longInt)atol(strReaction.c_str()),
                                         (acs_longInt)atol(strCounter.c_str()), (acs_double)atof(strKcond.c_str()), (acs_double)atof(strKcleav.c_str()),
                                         (acs_double)atof(strKcpx.c_str()), (acs_int)atoi(strCpxTarget.c_str())));

    }

    myfile.close();
    inserSubListInSpecies();
    if(debugLevel == FINDERRORDURINGRUNTIME) cout << "environment::createInitialCatalysisLayerFromSpecificFileSTD end" << endl;
    return true;
}//eof createInitialCatalysisLayerFromSpecificFileSTD


/**
 Create the diffusion constant renforcement according to the species length
 @version 1.0
 @param tmpDiffEnh diffusion enhancement from parameters
 @param MTRand& tmp_RndDoubleGen random number generator
 @param tmpNewSpeciesLength Lenght of the species
 */
acs_double environment::createDiffusionRenforcement(acs_double tmpDiffEnh, acs_int tmpNewSpeciesLength)
{
	if(debugLevel == FINDERRORDURINGRUNTIME) cout << "\t\t\t\tenvironment::createDiffusionRenforcement start" << endl;

	// THEORETICALLY DIFFUSION GOES WITH THE ROOT OF THE LENGTH, SO DIFFUSION KINETIC CONSTANT
	// IS EQUAL TO L^(-TMPDIFFENH)
	acs_double tmpDiffusionRenforcement = pow((double)tmpNewSpeciesLength, -(double)tmpDiffEnh);
	
	if(debugLevel == FINDERRORDURINGRUNTIME) cout << "\t\t\t\tenvironment::createDiffusionRenforcement end" << endl;

	
	return tmpDiffusionRenforcement;
}

/**
 Create the precipitation constant renforcement according to the species length
 @version 2.5.1
 @param tmpPreEnh precipitation enhancement from parameters
 @param MTRand& tmp_RndDoubleGen random number generator
 @param tmpNewSpeciesLength Lenght of the species
 */
bool environment::setSolubility(acs_int tmpNewSpeciesLength, MTRand& tmpRndDoubleGen)
{
	if(debugLevel == FINDERRORDURINGRUNTIME) cout << "\t\t\t\tenvironment::setSolubility start" << endl;

	bool tmpSolubility = true;
//	acs_double tmpDice;
//	acs_double lengthToLogistic;
	
        if(solubility_threshold > 0)
	{	
//		tmpDice = tmpRndDoubleGen();
//		lengthToLogistic = 1 / (1 + pow(NEP, (double)solubility_threshold-tmpNewSpeciesLength));
//		if(tmpDice < lengthToLogistic)
//		{
//			tmpSolubility = false;
//		}
                // COMMENT THIS LINE IF YOU WANT TO REACTIVATE THE TANH FUNCTION
                if(tmpNewSpeciesLength >= solubility_threshold)
                    tmpSolubility = false;
	}
	
	if(debugLevel == FINDERRORDURINGRUNTIME) cout << "\t\t\t\tenvironment::setSolubility end" << endl;
	
	return tmpSolubility;
}

/**
 This functions returns the ID of the species whether this one is already present, otherwise the new ID is created as the 
 number of elements in the species vector
 @version 1.0
 @param tmpNewSequence New sequence to evaluate
 @param MTRand& tmp_RndDoubleGen random number generator
 @param tmpNewSpeciesLength Lenght of the species
 */
acs_longInt environment::returnPosSpeciesAlreadyPresent(string tmpNewSequence)
{
	if(debugLevel == FINDERRORDURINGRUNTIME) cout << "\t\t\t\t\tenvironment::returnPosSpeciesAlreadyPresent start" << endl;

	acs_longInt tmpPosition = (acs_longInt)allSpecies.size();
	for(acs_longInt i = 0; i < getNumberOfTheoreticalSpecies(); i++)
	{
		if(allSpecies.at(i).getSequence() == tmpNewSequence)
		{
			if(allSpecies.at(i).getComplexCutPnt() == 0)
			{
				tmpPosition = i;
				break;
			}
		}
	}	
	
	if(debugLevel == FINDERRORDURINGRUNTIME) cout << "\t\t\t\t\tenvironment::returnPosSpeciesAlreadyPresent end" << endl;

	return tmpPosition;
}

/**
 Return the reaction ID if the reaction is already present, otherwise it returns the new reaction ID
 @version 1.0
 @param acs_int tmpReactionType rection type (cleavage or condensation)
 @param acs_longInt tmpIds_I species I ID
 @param acs_longInt tmpIds_II species II ID
 @param acs_longInt tmpIds_III species III ID
 */
acs_longInt environment::returnPosReactionAlreadyPresent(acs_int tmpReactionType, acs_longInt tmpIds_I, acs_longInt tmpIds_II, acs_longInt tmpIds_III)
{
	if(debugLevel == FINDERRORDURINGRUNTIME) cout << "\t\t\t\t\tenvironment::returnPosReactionAlreadyPresent start" << endl;

	acs_longInt reactionPos = (acs_longInt)allReactions.size();
	for(acs_longInt i = 0; i < getNumberOfReactions(); i++) //FOR EACH REACTION
	{
		if(allReactions.at(i).getType() == tmpReactionType) //IF THE TYPE IS THE SAME
		{
			if(allReactions.at(i).getSpecies_I() == tmpIds_I) //IF THE FIRST SPECIES IS THE SAME
			{
				// IF THE SECOND AND THE THIRD SPECIES ARE THE SAME, THE POSITION DOESN'T MATTER 
				if( (allReactions.at(i).getSpecies_II() == tmpIds_II && allReactions.at(i).getSpecies_III() == tmpIds_III) ||
					(allReactions.at(i).getSpecies_II() == tmpIds_III && allReactions.at(i).getSpecies_III() == tmpIds_II)) 
				{
					reactionPos = i;
					break;
				}
			}
		}	
	}
	
	if(debugLevel == FINDERRORDURINGRUNTIME) cout << "\t\t\t\t\tenvironment::returnPosReactionAlreadyPresent end" << endl;
	
	return reactionPos;
}

/**
 If the reaction is not new this function checks if the reactions has been already catalysed by this species
 @version 1.0
 @param acs_int tmpSPeciesID catalyst ID
 @param tmpIdReaction reaction ID
 */
bool environment::checkIfTheReactionIsAlreadyCatalyzedByThisSpecies(acs_longInt tmpSPeciesID, acs_longInt tmpIdReaction)
{
	if(debugLevel == FINDERRORDURINGRUNTIME) cout << "\t\t\t\t\tenvironment::checkIfTheReactionIsAlreadyCatalyzedByThisSpecies start" << endl;

	bool tmpAlreadyCatFlag = false;
	for(acs_longInt i = 0; i < getNumberOfCatalysis(); i++)
	{
		if(allCatalysis.at(i).getCat() == tmpSPeciesID)
		{
			if(allCatalysis.at(i).getReactionID() == tmpIdReaction)
			{
				tmpAlreadyCatFlag = true;
				break;
			}
		}
	}
	
	if(debugLevel == FINDERRORDURINGRUNTIME) cout << "\t\t\t\t\tenvironment::checkIfTheReactionIsAlreadyCatalyzedByThisSpecies end" << endl;
	
	return tmpAlreadyCatFlag;

}

/**
 This function return false if the tmpCandidateProduct is a catalyst of tmpCat
 @version 1.0 - last update 2009/10/08 - build 009
 
 */
bool environment::checkIfOnlyMutualCatalysis(acs_int tmpCat, acs_int tmpCandidateProduct)
{
	if(debugLevel == FINDERRORDURINGRUNTIME) cout << "\t\t\t\t\tenvironment::checkIfOnlyMutualCatalysis start" << endl;
	bool tmpNoMutualCatalysis = true;
	for(acs_longInt i = 0; i < getNumberOfCatalysis(); i++)
	{
		// Select reaction for the candidate product
		if(allCatalysis.at(i).getCat() == tmpCandidateProduct)
		{
			if(allReactions.at(allCatalysis.at(i).getReactionID()).getType() == 0) // Condensation
			{
				if(allReactions.at(allCatalysis.at(i).getReactionID()).getSpecies_I() == tmpCat)
				{
					tmpNoMutualCatalysis = false;
					break;
				}
			}else{	// Cleavage
				if(allReactions.at(allCatalysis.at(i).getReactionID()).getSpecies_II() == tmpCat || 
				   allReactions.at(allCatalysis.at(i).getReactionID()).getSpecies_III() == tmpCat)
				{
					tmpNoMutualCatalysis = false;
					break;
				}				
			}
		}
	}
	

	if(debugLevel == FINDERRORDURINGRUNTIME) cout << "\t\t\t\t\tenvironment::checkIfOnlyMutualCatalysis end" << endl;
	return tmpNoMutualCatalysis;
}


/**
 This function check the coherence of the data structures. IF IT FAULTS THE SYSTEM RETURNS ERROR!!!
 IT IS VERY RECOMMENDED RUN THIS FUNCTION ALWAYS AFTER THE INITIALIZATION
 @version 1.1
 @date 2011-04-13
 */
bool environment::structureCoherenceCheckUp()
{
	if(debugLevel == FINDERRORDURINGRUNTIME) cout << "environment::structureCoherenceCheckUp start" << endl;
	
	bool flagControl = true;
	
	if(debugLevel >= RUNNING_VERSION)
	{
		cout << endl;
		cout << "|- STRUCTURE COHERENCE CONTROL |" << endl;
	
	// ===================================
	// SPECIES STRUCTURE COHERENCE CONTROL
	// ===================================
	
		cout << "\t|- SPECIES STRUCTURE" << endl;
		cout << "\t\t|- Searching for Incoherences, double entries and incorrect ID progression... "<< endl;
	}
	if(debugLevel >= RUNNING_VERSION)
	{
		acs_int tmpTotTheorSpecies = getNumberOfTheoreticalSpecies();
		for(acs_longInt i = 0; i < tmpTotTheorSpecies; i++)
		{
			if(i%100 == 0) cout << "\t\t|- Species Percentage Controlled  " << (double)i/tmpTotTheorSpecies*100 << " %" << endl;
			for(acs_longInt j= 0; j < tmpTotTheorSpecies; j++)
			{
				if(i != j) // A Species is definitely equal to itself
				{
					if(allSpecies.at(i).getSequence() == allSpecies.at(j).getSequence())
					{ // IF THE SEQUENCE IS THE SAME THE CUTTING POINT HAS TO BE CHECKED
						if(debugLevel >= SMALL_DEBUG)
						{
							cout << endl << "\t\t\t |- The sequence of " <<  allSpecies.at(i).getSequence() << " [ID" << i << "] is the same of "
								 <<   allSpecies.at(i).getSequence() << " [ID" << j << "]... check CUTTING POINT";
						}
					   if(allSpecies.at(i).getComplexCutPnt() == allSpecies.at(j).getComplexCutPnt())
					   {
						   flagControl = false;
						   cout << endl << "\t\t\t|- ERROR: Species " << allSpecies.at(i).getID() << " is a duplicate of species "
								<< allSpecies.at(j).getID() << ": The two cutting points are NOT different: "
								<< allSpecies.at(i).getComplexCutPnt()  << " and " << allSpecies.at(j).getComplexCutPnt()  << endl;
						   break;
					   }else{
						   if(debugLevel >= SMALL_DEBUG)
						   {
							   cout << endl << "\t\t\t|- NO PROBLEN: The two cutting points are different: " << allSpecies.at(i).getComplexCutPnt()
									<< " and " << allSpecies.at(j).getComplexCutPnt()  << endl;
						   }
					   }
					}
				}
			}
			if(debugLevel >= RUNNING_VERSION)
			{
				if(i == getNumberOfTheoreticalSpecies())
				{
					if(flagControl == true)
						cout << endl << "\t\t\t|- No Double entries... OK!" << endl;
				}
			}

			if(debugLevel >= SMALL_DEBUG)
			{
				cout << endl << "ID: " << allSpecies.at(i).getID() << "\t"
				<< "SEQ: " << allSpecies.at(i).getSequence() << "\t"
				<< "AMOUNT: " << allSpecies.at(i).getAmount() << "\t"
				<< "DIFF: " << allSpecies.at(i).getDiffusionEnh() << "\t"
				<< "PREC: " << allSpecies.at(i).getSolubility() << "\t"
				<< "COM DEG: " << allSpecies.at(i).getComplexDegEnh() << "\t"
				<< "CUT PNT: " << allSpecies.at(i).getComplexCutPnt() << "\t"
				<< "EVAL: " << allSpecies.at(i).getEvaluated();
			}
					//SINGLE CHECKS, with the following TRY all allSpecies.at(i) are checked
					try{
						if(allSpecies.at(i).getID() != i) //CHECK ID PROGRESSION
						{
							flagControl = false;
							cout << endl << "\t\t\tERROR: SPECIES ID progression in not right: ID " << allSpecies.at(i).getID() << " could be equal to " << i << endl;
						}
					}catch(exception&e){
						cout << "structureCoherenceCheckUp allSpecies.at(i)…" << endl;
						cout << "Vectorsize " << allSpecies.size() << " - position " << i << endl;
						cerr << "exceptioncaught:" << e.what() << endl;
						ExitWithError("structureCoherenceCheckUp","exceptionerrorthrown");
				   }

			if(allSpecies.at(i).getComplexCutPnt() > 0 && allSpecies.at(i).getEvaluated() > 0) // COMPLEX CANNOT BE EVALUATED
			{
				flagControl = false;
				cout << endl << "Species " << allSpecies.at(i).getID() << " is a complex and cannot be evaluated!!!" << endl;
			}
			// IF THERE ARE MOLECULE BELONGING TO THE SPECIES, THE SPECIES HAD TO BE EVALUATED
			if(allSpecies.at(i).getAmount() > 0 &&  allSpecies.at(i).getComplexCutPnt() == 0 && allSpecies.at(i).getEvaluated() == 0)
			{
				flagControl = false;
				cout << endl << "Species " << allSpecies.at(i).getSequence() << "(ID" << i << ") has "
					 <<  allSpecies.at(i).getAmount() << " and it wasn't evaluated!!!" << endl;
			}
			// A COMPLEX HAS TO BE EOLUBILITY SET TO 1
			if(allSpecies.at(i).getComplexCutPnt() > 0 && allSpecies.at(i).getSolubility() == PRECIPITATED)
			{
				flagControl = false;
				cout << endl << "Complex " << allSpecies.at(i).getSequence() << "(ID" << i << ") has to be solubility set to 1!!!" << endl;
			}
		}
		
		if(debugLevel >= RUNNING_VERSION)
		{
			//Previous Test
			cout << endl << "\t\t\t|- ID Progression... OK!" << endl;
			cout << "\t\t\t|- Complexes not evaluated... OK!" << endl;
			cout << "\t\t\t|- Living Species Evaluated... OK!" << endl;
			cout << "\t\t\t|- Complex Solubility... OK!" << endl;
		}
	}
	
	// ===================================
	// IF THE SYSTEM IS OPEN
	// ===================================
	
    if(influx_rate > 0)
	{
            if(flagControl) // IF THE PREVIOUS CONTROL WAS OK REACTION CONTROL WILL BE PERFORMED
            {
                if(debugLevel == RUNNING_VERSION)
                {
                    //Check reactions structures
                    cout << "\t|- INFLUX LAYER STRUCTURE" << endl;
                    cout << "\t\t|- Searching for firing disk lacking coupling with the influx species and incorrect cumulative probability... ";
                }
                acs_int tempPresentSpeciesTOT = 0;
                acs_int tempInfluxID = 0;

                for(vector<acs_int>::iterator tmpAllNutrientsIter = nutrientsForInflux.begin(); tmpAllNutrientsIter != nutrientsForInflux.end(); tmpAllNutrientsIter++)
                {
                    for(vector<species>::iterator tmpAllSpeciesIter = allSpecies.begin(); tmpAllSpeciesIter != allSpecies.end(); tmpAllSpeciesIter++)
                    {
                        if(*tmpAllNutrientsIter == tmpAllSpeciesIter->getID())
                        {
                                tempPresentSpeciesTOT += 1;
                                break;
                        }
                    }
                    tempInfluxID += 1;
                }

                if(tempInfluxID != tempPresentSpeciesTOT)
                {
                    flagControl = false;
                    cout << endl << "\t\t\tERROR: Not all the influx species are prenset in the firing fisk";
                }else{
                    if(debugLevel == RUNNING_VERSION)
                    {
                            cout << endl << "\t\t\t|- Influx Species are already present in the firing disk... OK!" << endl;

                    }
                }
            }
	}
	
	// ===================================
	// IF THERE ARE REACTIONS
	// ===================================

	if(debugLevel >= RUNNING_VERSION)
	{
		if((acs_longInt)allReactions.size() > 0)
		{
				if(debugLevel == RUNNING_VERSION)
				{
					//Check reactions structures
					cout << "\t|- REACTIONS LAYER STRUCTURE" << endl;
					cout << "\t\t|- Searching for duplicate, incorrect reactions products and ID progression... " << endl;
				}
	
				if(flagControl) // IF THE PREVIOUS CONTROL WAS OK REACTION CONTROL WILL BE PERFORMED
				{
					string tmpStr1; // TEMPORARY TEST STRING 1
					string tmpStr2; // TEMPORARY TEST STRING 2
					string tmpStr3; // TEMPORARY TEST STRING 3
					string tmpStrtot; // TEMPORARY CONCATENATION OF STR2 AND STR3
					acs_double tmpClConRatio = 0;

					// FOR ALL REACTIONS
					acs_int tmpNrct = getNumberOfReactions();
					for(acs_longInt i = 0; i < tmpNrct - 1; i++)
					{
						if(i%100 == 0) cout << "\t\t|- Reaction Percentage Controlled " << (double)i/tmpNrct*100 << " %" << endl;
						if (allReactions.at(i).getType() == 1) {
								tmpClConRatio++;
						}

						if(debugLevel >= MEDIUM_DEBUG)
						{
							cout << endl << "ID: " << allReactions.at(i).getID() <<
											"\tType: " << allReactions.at(i).getType() <<
											"\tMOL 1: " << allReactions.at(i).getSpecies_I() <<
											"\tMOL 2: " << allReactions.at(i).getSpecies_II() <<
											"\tMOL 3: " << allReactions.at(i).getSpecies_III() <<
											"\tOccurence: " << allReactions.at(i).getEvents() <<
											"\tEquilibrium constant: " << allReactions.at(i).getEnergyType();
						}

						// CHECK ID PROGRESSION
						if(allReactions.at(i).getID() != i)
						{
							flagControl = false;
							cout << endl << "\t\t\tERROR: REACTION ID progression in not right: ID " << allReactions.at(i).getID() <<
							" could be equal to " << i << endl;
							//break;
						}

						// CHECK THAT SUBSTRATES ARE NOT INVOLVED IN REACTIONS
						try{
							if(allSpecies.at(allReactions.at(i).getSpecies_I()).getComplexCutPnt() > 0)
							{
									flagControl = false;
									cout << endl << "\t\t\tERROR: Reaction " << i << ", Species " << allReactions.at(i).getSpecies_I() << " ia a complex, then it can not be affected by any reaction!!!" << endl;
							}
						}catch(exception&e){
							cout << "structureCoherenceCheckUp allSpecies.at(allReactions.at(i).getSpecies_I())…" << endl;
							cout << "Vectorsize " << allSpecies.size() <<" - position "<< i << endl;
							cerr << "exceptioncaught:"<<e.what()<<endl;
							ExitWithError("structureCoherenceCheckUp","exceptionerrorthrown");
					   }

					   try{
							if(allSpecies.at(allReactions.at(i).getSpecies_II()).getComplexCutPnt() > 0)
							{
									flagControl = false;
									cout << endl << "\t\t\tERROR: Reaction " << i << ", Species " << allReactions.at(i).getSpecies_II() << " ia a complex, then it can not be affected by any reaction!!!" << endl;
							}
						}catch(exception&e){
							cout << "if(allSpecies.at(allReactions.at(i).getSpecies_II()).getComplexCutPnt() > 0)" << endl;
							cout << "Vectorsize "<<allSpecies.size()<<" - position " << i << endl;
							cerr << "exceptioncaught:"<<e.what()<<endl;
							ExitWithError("structureCoherenceCheckUp","exceptionerrorthrown");
					   }

					   try{
							if(allSpecies.at(allReactions.at(i).getSpecies_III()).getComplexCutPnt() > 0)
							{
									flagControl = false;
									cout << endl << "\t\t\tERROR: Reaction " << i << ", Species " << allReactions.at(i).getSpecies_III() << " ia a complex, then it can not be affected by any reaction!!!" << endl;
							}
						}catch(exception&e){
							cout << "if(allSpecies.at(allReactions.at(i).getSpecies_III()).getComplexCutPnt() > 0)" << endl;
							cout << "Vectorsize "<<allSpecies.size()<<" - position " << i << endl;
							cerr << "exceptioncaught:"<<e.what()<<endl;
							ExitWithError("structureCoherenceCheckUp","exceptionerrorthrown");
					   }

						// CHECK COMPOSITION COHERENCE
						tmpStr1 = allSpecies.at(allReactions.at(i).getSpecies_I()).getSequence();
						tmpStr2 = allSpecies.at(allReactions.at(i).getSpecies_II()).getSequence();
						tmpStr3 = allSpecies.at(allReactions.at(i).getSpecies_III()).getSequence();
						tmpStrtot = tmpStr2 + tmpStr3;

						if(tmpStr1 != tmpStrtot)
						{
								flagControl = false;
								if(debugLevel >= RUNNING_VERSION)
										cout << endl << "\t\t\tERROR in reaction " << i << " - Type: " << allReactions.at(i).getType() << endl;
										cout << "\t\t\t\t: " << tmpStr1 << "[" << allReactions.at(i).getSpecies_I() << "] cannot be product or substrate of " << tmpStr2 << "[" << allReactions.at(i).getSpecies_II()
										<< "] + " << tmpStr3 << "[" << allReactions.at(i).getSpecies_III() << "]";
						}

						// SEARCHING FOR DUPLICATE REACTIONS
						for(acs_longInt j = 0; j < getNumberOfReactions(); j++)
						{
							if( i != j)// A single reaction is definitely equal to itself
							{
								try{
									if(allReactions.at(i).getType() == allReactions.at(j).getType())
									{
										if(allReactions.at(i).getSpecies_I() == allReactions.at(j).getSpecies_I())
										{
											if((allReactions.at(i).getSpecies_II() == allReactions.at(j).getSpecies_II() &&
												allReactions.at(i).getSpecies_III() == allReactions.at(j).getSpecies_III()))
											{
												flagControl = false;
												cout << endl << "\t\t\t Reaction " << allReactions.at(i).getID() << " " << allReactions.at(i).getType()
													 << " " << allReactions.at(i).getSpecies_I() << " " <<  allReactions.at(i).getSpecies_II()
													 << " " << allReactions.at(i).getSpecies_III() << endl;
												cout << endl << "\t\t\t Reaction " << allReactions.at(j).getID() << " " << allReactions.at(j).getType()
													 << " " << allReactions.at(j).getSpecies_I() << " " <<  allReactions.at(j).getSpecies_II()
													 << " " << allReactions.at(j).getSpecies_III() << endl;
												cout << endl << "\t\t\tERROR: Reaction " << allReactions.at(i).getID() << " is a duplicate of reaction " << allReactions.at(j).getID();
												break;
											}
										} // end if(allReactions.at(i).getSpecies_I() == allReactions.at(j).getSpecies_I())
									} // end if(allReactions.at(i).getType() == allReactions.at(j).getType())
								}catch(exception&e){
									cout << "allReactions.at(i).get..." << endl;
									cout << "Vectorsize "<<allReactions.size()<<" - position " << i << endl;
									cerr << "exceptioncaught:"<<e.what()<<endl;
									ExitWithError("structureCoherenceCheckUp","exceptionerrorthrown");
							   }
							} // end if( i != j)
						} // end for(acs_longInt j = 0; j < getNumberOfReactions(); j++)
					} // end for(acs_longInt i = 0; i < getNumberOfReactions() - 1; i++)
					cout << endl << "\t\t\t|- Ratio between cleavage and condensation reactions:" << tmpClConRatio / getNumberOfReactions() * 100;

				} // end if(flagControl)
		} // if((acs_longInt)allReactions.size() > 0)
	} // if(debugLevel >= MEDIUM_DEBUG)

	if(debugLevel >= RUNNING_VERSION)
	{
		// ===================================
		// IF THERE ARE CATALYSIS
		// ===================================

		if(flagControl)
		{
			if(debugLevel >= RUNNING_VERSION)
			{
				// PREVIOUS CHECK
				cout << endl
				<< "\t\t\t|- Id Progression... OK!" << endl
				<< "\t\t\t|- No Double reactions... OK!" << endl
				<< "\t\t\t|- Substrates are not involved in reactions... OK!" << endl
				<< "\t\t\t|- Reaction coherence... OK!" << endl;
				if (energy == 1) {
						cout << "\t\t\t|- ENERGY TYPE control... OK!" << endl;
				}
				// CATALYSIS STRUCTURE CHECK
				cout << "\t|- CATALYSIS LAYER STRUCTURE" << endl;
				cout << "\t\t|- Searching for duplicate, incorrect ID progression... " << endl;
			}
			acs_int tmpCatNum = getNumberOfCatalysis();
			if(tmpCatNum > 0)
				{
				for(acs_longInt i = 0; i < tmpCatNum - 1; i++)
				{
					if(i%100 == 0) cout << "\t\t|- Catalysis Percentage Controlled " << (double)i/tmpCatNum*100 << " %" << endl;
					// NO COMPLEX CATALYST CONTROL
					if(allSpecies.at(allCatalysis.at(i).getCat()).getComplexCutPnt() > 0)
					{
						//A COMPLEX CANNOT CATALYZE ANY REACTION
						flagControl = false;
						cout << "Catalysis number " << i << ": A complex cannot catalyse any reaction, catalyst " << allCatalysis.at(i).getCat() << " is a complex!!! " << endl;
						break;
					}
					// NO DOUBLE ENTRIES CONTROL
					for(acs_longInt j= 0; j < getNumberOfCatalysis(); j++)
					{
						if(i != j) // A Species is definitely equal to itself
						{
							if(allCatalysis.at(i).getCat() == allCatalysis.at(j).getCat())
							{
								if(allCatalysis.at(i).getReactionID() == allCatalysis.at(j).getReactionID())
								{
									flagControl = false;
									cout << endl << "\t\t\tERROR: Catalysis " << allCatalysis.at(i).getCatId() << " is a duplicate of catalysis " << allCatalysis.at(j).getCatId()  << endl;
									break;
								}
							}
						}
					}
					// PROGRESSION CONTROL
					if(allCatalysis.at(i).getCatId() != i)
					{
						flagControl = false;
						cout << endl << "\t\t\tERROR: CATALYSIS ID progression in not right: ID " << allCatalysis.at(i).getCatId() <<
						" could be equal to " << i << endl;
						//break;
					}

					// NOT INVERSE REACTIONS CONTROL
					if(main_rev_rct_allowed == false)
					{
						if(allReactions.at(allCatalysis.at(i).getReactionID()).getType() == CLEAVAGE)
						{ // CLEAVAGE
								if(!notInverseReactionAlreadyCatalyzed(CLEAVAGE,
																	   allReactions.at(allCatalysis.at(i).getReactionID()).getSpecies_I(),
																	   allReactions.at(allCatalysis.at(i).getReactionID()).getSpecies_II()))
								{
										flagControl = false;
										cout << endl << "\t\t\tERROR: Rct Type: " <<  CLEAVAGE  <<
										" - Catalisys: " << i << " - Catalyst: " <<  allCatalysis.at(i).getCat() <<
										" - S1: " << allReactions.at(allCatalysis.at(i).getReactionID()).getSpecies_I() <<
										" - S2: " << allReactions.at(allCatalysis.at(i).getReactionID()).getSpecies_II() <<
										" - Reaction " << allCatalysis.at(i).getReactionID() << endl << "\t\t\t\t" <<
										": Both forward and backward reaction are catalyzed!!!" << endl <<
										"\t\t\t\t" << "-----------------------------------------------------" << endl;
										//break;
								}
							// NO SELF-CATALYSIS CONTROL
		//					if(allCatalysis.at(i).getCat() < (acs_longInt)firingDisk.size())
		//					{
		//						if((allCatalysis.at(i).getCat() == allReactions.at(allCatalysis.at(i).getReactionID()).getSpecies_II()) ||
		//						   (allCatalysis.at(i).getCat() == allReactions.at(allCatalysis.at(i).getReactionID()).getSpecies_III()))
		//						{
		//							cout << "WARNING*** CLEAVAGE -> Species " << allCatalysis.at(i).getCatId() << " is a firing disk species and it cannot catalize the formation of itself: cat ID "
		//							     << i << " - Rct ID " << allCatalysis.at(i).getReactionID() << endl;
		//							//flagControl = false;
		//						}
		//					}
						}else{ // CONDENSATION
								if(!notInverseReactionAlreadyCatalyzed(CONDENSATION,
													   allReactions.at(allCatalysis.at(i).getReactionID()).getSpecies_I(),
													   allReactions.at(allCatalysis.at(i).getReactionID()).getSpecies_II()))
								{
										flagControl = false;
										cout << endl << "\t\t\tERROR: Rct Type: " << CONDENSATION <<
										" - Catalisys: " << i << " - Catalyst: " <<  allCatalysis.at(i).getCat() <<
										" - S1: " << allReactions.at(allCatalysis.at(i).getReactionID()).getSpecies_I() <<
										" - S2: " << allReactions.at(allCatalysis.at(i).getReactionID()).getSpecies_II() <<
										" - Reaction " << allCatalysis.at(i).getReactionID() << endl << "\t\t\t\t" <<
										": Both forward and backward reaction are catalyzed!!!" << endl <<
										"\t\t\t\t" << "-----------------------------------------------------" << endl;
										//break;
								}
							// NO SELF-CATALYSIS CONTROL
		//					if(allCatalysis.at(i).getCat() < (acs_longInt)firingDisk.size())
		//					{
		//						if((allCatalysis.at(i).getCat() == allReactions.at(allCatalysis.at(i).getReactionID()).getSpecies_I()))
		//						{
		//							cout << "WARNING*** CONDENSATION -> Species " << allCatalysis.at(i).getCatId() << " is a firing disk species and it cannot catalize the formation of itself: cat ID "
		//							<< i << " - Rct ID " << allCatalysis.at(i).getReactionID() << endl;
		//							//flagControl = false;
		//						}
		//					} // end if(allCatalysis.at(i).getCatId() < (acs_longInt)firingDisk.size())
						} // end if(allReactions.at(allCatalysis.at(i).getReactionID()).getType() == CLEAVAGE)
					} // end if(main_rev_rct_allowed == false)
				} // end for(acs_longInt i = 0; i < getNumberOfCatalysis() - 1; i++)
			}else{ // end if(tmpCatNum > 0) If there are catalysis
				cout << "\t\t\t	CATALYSED REACTIONS NOT PRESENT!!! " << CONDENSATION << endl;
			}
		} // end if(flagControl)
	}
	
	if(flagControl)
	{
		if(debugLevel >= RUNNING_VERSION)
		{
			cout << endl
				<< "\t\t\t|- No doubles in catalysis structures... OK!" << endl
				<< "\t\t\t|- Id Progression... OK!" << endl
				<< "\t\t\t|- No inverse reactions are catalyzed by the same species... OK!" << endl
				<< "\t\t\t|- Reactions not catalyzed by complexes... OK!" << endl
				<< "\t\t\t|- No Autocatalysis into the firing disk... OK!" << endl;
		}
	}
	
	if(debugLevel == FINDERRORDURINGRUNTIME) cout << "environment::structureCoherenceCheckUp end" << endl;
	
	return flagControl;

}

/**
 Insert substrate reactions information in species if complex
 @version 1.0
 @date 20130917
 @author Alessandro Filisetti
 */
void environment::inserSubListInSpecies()
{
	if(debugLevel == FINDERRORDURINGRUNTIME) cout << "environment::inserSubListInSpecies start" << endl;
	// For Each Species
	acs_longInt tmpCat;
	acs_longInt tmpSub;
	acs_longInt tmpRct;
	acs_int tmpCpxTgt;

	if(debugLevel == RUNNING_VERSION)
	{
		cout << "\t|- Complex Condensation Parameters Introduction...";
	}

    for(vector<species>::iterator tmpSpeciesIterator = allSpecies.begin(); tmpSpeciesIterator != allSpecies.end(); tmpSpeciesIterator++)
    {
    	// If the species is a complex
    	if(tmpSpeciesIterator->getComplexCutPnt() > 0)
    	{
    		tmpCat = tmpSpeciesIterator->getCatalyst_ID();
    		tmpSub = tmpSpeciesIterator->getSubstrate_ID();
    		// For each catalysis
    		for(vector<catalysis>::iterator tmpCatalysisIterator = allCatalysis.begin(); tmpCatalysisIterator != allCatalysis.end(); tmpCatalysisIterator++)
    		{
    			// Searching for Catalysed reactions
    			if(tmpCatalysisIterator->getCat() == tmpCat)
    			{
    				tmpRct = tmpCatalysisIterator->getReactionID();
    				tmpCpxTgt = tmpCatalysisIterator->getCpxTarget();
    				// Insert complex formation reaction
    				if((allReactions.at(tmpRct).getType()) == CONDENSATION || (reverseReactions == true))
    				{
    					if(debugLevel == SMALL_DEBUG)
    					{
    						cout << tmpSpeciesIterator->getID() << " " << tmpCatalysisIterator->getCatId() << endl;
    						cout << tmpSpeciesIterator->getCatalyst_ID() << " " << allReactions.at(tmpRct).getID() << endl;
    					}
						if(tmpCpxTgt == 1)
						{
							if(allReactions.at(tmpRct).getSpecies_II() == tmpSub)
								tmpSpeciesIterator->insertSecSub(allReactions.at(tmpRct).getSpecies_III(),tmpCatalysisIterator->getKass(),tmpCatalysisIterator->getCatId());
							//cout << allReactions.at(tmpRct).getSpecies_III() << " " << tmpCatalysisIterator->getKass() << endl;
						}else{
							if(allReactions.at(tmpRct).getSpecies_III() == tmpSub)
								tmpSpeciesIterator->insertSecSub(allReactions.at(tmpRct).getSpecies_II(),tmpCatalysisIterator->getKass(),tmpCatalysisIterator->getCatId());
							//cout << allReactions.at(tmpRct).getSpecies_II() << " " << tmpCatalysisIterator->getKass() << endl;
						}
    				}
    			}
    		}
    	}
    }

    //showSubListInSpecies();

    if(debugLevel == FINDERRORDURINGRUNTIME) cout << "environment::inserSubListInSpecies end" << endl;
}

/**
 * Function to show the substrates list of the complex species
 */
void environment::showSubListInSpecies()
{
	if(debugLevel == FINDERRORDURINGRUNTIME) cout << "environment::showSubListInSpecies start" << endl;
	if(debugLevel == RUNNING_VERSION)
	{
		cout << "|- Complex Condensation Parameters ..." << endl;
	}

	for(vector<species>::iterator tmpSpeciesIterator = allSpecies.begin(); tmpSpeciesIterator != allSpecies.end(); tmpSpeciesIterator++)
	{
		// If the species is a complex
		if(tmpSpeciesIterator->getComplexCutPnt() > 0)
		{
			if(tmpSpeciesIterator->getSecSubListSize() > 0)
			{
				cout << "\t|- Complex " << tmpSpeciesIterator->getID() << " - " << tmpSpeciesIterator->getSequence()
					 << " - " << allSpecies.at(tmpSpeciesIterator->getCatalyst_ID()).getSequence()
					 << " - " << allSpecies.at(tmpSpeciesIterator->getSubstrate_ID()).getSequence() << endl;
				cout << "\t\t|- Sec Subs ";
				for(acs_int i = 0; i < tmpSpeciesIterator->getSecSubListSize(); i++)
				{
					cout << tmpSpeciesIterator->getSecSubListID(i) << " " << tmpSpeciesIterator->getSec_k_SubListID(i) << " || ";
				}
				cout << endl;
			}
		}
	}
	if(debugLevel == FINDERRORDURINGRUNTIME) cout << "environment::showSubListInSpecies end" << endl;
}

/**
 * Function to show the Gillespie engagement of each species.
 * @version: 1.0
 * @date: 20131031
 * @author: Alessandro Filisetti
 */
void environment::showGillEngagementInSpecies()
{
	if(debugLevel == FINDERRORDURINGRUNTIME) cout << "environment::showGillEngagementInSpecies start" << endl;
	if(debugLevel == RUNNING_VERSION)
	{
		cout << "|- Gillespie species engagement ..." << endl;
	}

	for(vector<species>::iterator tmpSpeciesIterator = allSpecies.begin(); tmpSpeciesIterator != allSpecies.end(); tmpSpeciesIterator++)
	{
		tmpSpeciesIterator->showGillEngagement();
	}
	if(debugLevel == FINDERRORDURINGRUNTIME) cout << "environment::showGillEngagementInSpecies end" << endl;
}

/**
 Check if the reaction catalyze both one reaction and the inverted one
 @version 1.0
 */
bool environment::notInverseReactionAlreadyCatalyzed(acs_int tmpRct, acs_longInt tmpID_I, acs_longInt tmpID_II)
{
	if(debugLevel == FINDERRORDURINGRUNTIME) cout << "environment::notInverseReactionAlreadyCatalyzed start" << endl;
	
	// BOOLEAN VARIABLE CONTROL
	bool reactionFlag = true;
	
	int tmpRevRctType = CLEAVAGE;
	if(tmpRct == CLEAVAGE)
		tmpRevRctType = CONDENSATION;
	
	for(acs_int tmpMMRctId = 0; tmpMMRctId < (acs_longInt)allReactions.size(); tmpMMRctId++)
	{
            try{
                if(allReactions.at(tmpMMRctId).getType() == (unsigned)tmpRevRctType)
                {
                    //IF THE INVERTED REACTION IS ALREADY CATALYSED
                    if((allReactions.at(tmpMMRctId).getSpecies_I() == tmpID_I) &&
                       (allReactions.at(tmpMMRctId).getSpecies_II() == tmpID_II))
                    {
                        if(debugLevel >= SMALL_DEBUG)
                        {
                            cout << endl <<
                            "\t\t\t *** WARNING: The system is trying to create a reverse reaction of type " << tmpRct <<
                            " with S1 " << tmpID_I << " and S2 " << tmpID_II << endl <<
                            "\t\t\t\t but in the system is already present the reverse reaction with ID " << tmpMMRctId <<
                            " of type " << allReactions.at(tmpMMRctId).getType() <<
                            " with S1 " << allReactions.at(tmpMMRctId).getSpecies_I() <<
                            " and S2 " << allReactions.at(tmpMMRctId).getSpecies_II() << endl;
                        }
                        reactionFlag = false;
                        break;
                    }
                }
            }catch(exception&e){
                cout << "allReactions.at(tmpMMRctId).getType() == (unsigned)tmpRevRctType)" << endl;
                cout << "Vectorsize "<<allSpecies.size()<<" - position " << tmpMMRctId << endl;
                cerr << "exceptioncaught:"<<e.what()<<endl;
                ExitWithError("notInverseReactionAlreadyCatalyzed","exceptionerrorthrown");
           }
	}	
	
	if(debugLevel == FINDERRORDURINGRUNTIME) cout << "environment::notInverseReactionAlreadyCatalyzed end" << endl;

	return reactionFlag;
}

/**
  Return TRUE if there are sufficient molecules for the reaction. It is used for the reaction in which catalyst and substrate are the same
  molecules
  @version 1.0
  @date 2011.07.25
  @param tmpMI
  @param tmpMII
  @param tmpQI
  @param tmpQII
  */
bool environment::checkAvailability(acs_longInt tmpMI, acs_longInt tmpMII, acs_longInt tmpQI, acs_longInt tmpQII)
{
    bool flagVar = false;
    if(tmpMI == tmpMII)
    {
        if(debugLevel == HIGH_DEBUG)
        {
            cout << "\t\t|- Catalyst and Substrate are the same, ID#" << tmpMI << " " << tmpMII << endl;
            cout << "\t\t|- Amounts --> " << tmpQI << " " << tmpQII << endl;
        }
        if(tmpQI != tmpQII)
            ExitWithError("checkAvailability", "!!! The species are the same, ALSO THE QUANTITIES MUST BE THE SAME");
        if(tmpQI >= 2)
            flagVar = true;
    }else{
        flagVar = true;
    }


    return flagVar;
}

/**
 Perform all the gillespie algorithm procedure (OPTIMIZED VERSION)
 @version 1.0
 @date 2013.09.17
 @param tmpRndDoubleGen random numbers generator
 @param tmpTimeElapsed Computation time elapsed from the T0
 @param tmpActGEN actual generation
 @param tmpActSIM actual simulation
 @param tmpActSTEP actual step (reaction)
 @param tmpStoringPath path where results are stored
*/
bool environment::performOPTGillespieComputation(MTRand& tmpRndDoubleGen, clock_t& tmpTimeElapsed, acs_int tmpActGEN, acs_int tmpActSIM,
												 acs_int tmpActSTEP, string tmpStoringPath)
{
	if(debugLevel == FINDERRORDURINGRUNTIME) cout << "environment::performOPTGillespieComputation start" << endl;

	bool flagControl = true; //Check flag
	bool sameSpecies = false; // true if I'm computing the score of a reaction involving the same species as catalyst and substrate
	//acs_double temp_score; // Gillespie score
    gillespieTotalScore = 0; // Initialize (or reset) the overall Gillespie score
    gillespieNewSpeciesScore = 0; // Initialize (or reset) the sum of the Gillespie scores creating new species
    ratioBetweenNewGillTotGill = 0;
    reverseReactionsGillScore = 0; // the sum of the gillespie scores related to the reverse reactions
    ratioBetweenReverseAndTotalScore = 0; // Ratio between the reverse reaction scores and all the reaction scores
	acs_longInt temp_mol_I; // temp_mol_I II III and IV will contain the molecules IDs according to the different reactions
	acs_longInt temp_mol_II;
	acs_longInt temp_mol_III;
	acs_longInt temp_mol_IV;

	acs_longInt temp_substrateAmount; // These variables will contain the different amounts of the different species
	acs_longInt temp_cpxAmount;

	acs_longInt temp_catAmount_charged;
	acs_longInt temp_catAmount_NotCharged;
	acs_longInt temp_catAmount_TOT;

	acs_longInt temp_substrateAmount_charged;
	acs_longInt temp_substrateAmount_NotCharged;
	acs_longInt temp_substrateAmount_TOT;

	acs_int temp_rctType; // Reaction type (condensation, endo_condensation, complex_formation, endo_complex_formation, cleavage, endo_cleavage)
	acs_longInt temp_reactionID;
	acs_longInt temp_catalysisID;
	acs_longInt reaction_u; // This variable will contain the ID of the reaction selected amond all the possible one contained in the Gillespie algorithm
	acs_double tempTime; // Temporary Delta T
	acs_longInt tmpCatalyst_ID; //temporary catalyst produced degrading the complex to searching for the substrate in condensation reaction
	acs_longInt tmpSubstrate_ID; //temporary substrate produced degrading the complex to searching for the second substrate of the condensation reaction
	acs_double tmpDeltaT; // Time interval between two successive generation

	string nrgBooleanFunction = ""; // This string contains the energetic boolean function of the selected reaction

	acs_longInt tempIdCatalysis; // This variable containing the temporary ID of the catalysis

	// TRANSFORM KINETIC CONSTANTS IN GILLESPIE KINETIC CONSTANTS (THEY ARE RESCALED ACCORDING TO THE VOLUME)

	/* ------------------------------------------\
	// GILLESPIE STRUCTURE CREATION				 ||==----
	// ------------------------------------------*/

	// FOR EACH SINGLE SPECIES ALL POSSIBLE REACTIONS INVOLVING IT IS COMPUTED.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-

	// Initialize Time for the Gillespie Algorithm
    gillespiePartialTimer = clock();

    acs_longInt complexID;
    if(tmpActSTEP > 0) // The reaction is now created just in the first reaction (once new gill computation is finished condition will be == 1)
    {
    	if(debugLevel == SMALL_DEBUG) cout << "\t\t|- GILLESPIE STRUCTURE CREATION..." << endl;

    	// FOR EACH SPECIES
		for(vector<species>::iterator speciesIter = allSpecies.begin(); speciesIter != allSpecies.end(); speciesIter++)
		{
			complexID = speciesIter->getID();

			sameSpecies = false; // This flag will be true is a second order reaction will involve the same species

		// UPDATE SPECIES AGE .-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-

			if((speciesIter->getAmount() > 0) && (speciesIter->getComplexCutPnt() == 0))
			{
				// tmpActSTEP is the simulation step
				if(tmpActSTEP > 2)
				{
					speciesIter->setNewAge(reactionsTime.at((reactionsTime.size() - 1)) - reactionsTime.at((reactionsTime.size() - 2)));
				}
			}else{
				speciesIter->setNewAge(0); // IF THE SPECIES IS DEAD THE AGE OF LAST DEATH REMAINS
			}

			if((speciesIter->getAmount() > 0) && (speciesIter->getSolubility() == SOLUBLE) && (speciesIter->getComplexCutPnt() > 0)) // If there are some molecule belonging to this species
			{
				tmpCatalyst_ID = speciesIter->getCatalyst_ID();
				tmpSubstrate_ID = speciesIter->getSubstrate_ID();

				if(speciesIter->getSecSubListSize() > 0)
				{
					// If so, searching for the catalysed reactions involving the substrate
					// Check whether the temporary catalyst catalyse a condensation reaction containing the temporary substrate as a first or second molecule
					// WHETHER IT IS A CONDENSATION REACTION

					for(acs_int listCondSecStep = 0; listCondSecStep < speciesIter->getSecSubListSize(); listCondSecStep++)
					{
						// CONDENSATION SECOND STEP
						// Assign species to the temp variables
						temp_mol_I = speciesIter->getCatalyst_ID(); //catalyst
						temp_mol_II = speciesIter->getSecSubListID(listCondSecStep); // second substrate
						temp_catalysisID = speciesIter->getCatalysisIfCpxID(listCondSecStep); // Catalysis ID
						temp_mol_III = allReactions.at(allCatalysis.at(temp_catalysisID).getReactionID()).getSpecies_I(); // Product
						temp_mol_IV = complexID; // Complex
						temp_reactionID = allCatalysis.at(temp_catalysisID).getReactionID(); // reaction ID
						temp_rctType = CONDENSATION; // reaction type

						if(debugLevel == COMPLEXSTUFF)
						{
							if(temp_mol_I == 100 && temp_mol_II == 2)
							{
								cout << "**************************" << endl
									<< speciesIter->getCatalysisIfCpxID(listCondSecStep) << " " << temp_catalysisID << " " << temp_reactionID << " "
									<< allReactions.at(temp_reactionID).getSpecies_I()
									<< " " << allReactions.at(temp_reactionID).getSpecies_II()
									<< " " << allReactions.at(temp_reactionID).getSpecies_III()
									<< endl << "**************************" << endl
									<< endl;
							}
						}

						try{
							if(nrgBoolFlag == ENERGYBASED)
							{
								// COMPLEX CHARGED - SECOND SUBSTRATE NOT CHARGED
								if((nrgBooleanFunction[1] == TRUENRG) || (nrgBooleanFunction[3] == TRUENRG) || (nrgBooleanFunction[5] == TRUENRG))
								{
									// Compute the total amount of complexes and second substrates
									temp_cpxAmount = speciesIter->getChargeMols();
									temp_substrateAmount = allSpecies.at(temp_mol_II).getNOTchargeMols();
									performSingleGilleSpieIntroduction(temp_cpxAmount, temp_substrateAmount, complexID, temp_mol_II, temp_catalysisID, ENDO_CONDENSATION,
																	   temp_mol_I, temp_mol_II, temp_mol_III, temp_mol_IV, COMPLEXLOAD, temp_reactionID,false);
								}

								// COMPLEX NOT CHARGED - SECOND SUBSTRATE CHARGED
								if(nrgBooleanFunction[6] == TRUENRG)
								{
									// Compute the total amount of complexes and second substrates
									temp_cpxAmount = speciesIter->getNOTchargeMols();
									temp_substrateAmount = allSpecies.at(temp_mol_II).getChargeMols();
									performSingleGilleSpieIntroduction(temp_cpxAmount, temp_substrateAmount, complexID, temp_mol_II, temp_catalysisID, ENDO_CONDENSATION,
																	   temp_mol_I, temp_mol_II, temp_mol_III, temp_mol_IV, SUBSTRATELOAD, temp_reactionID,false);
								}

								// COMPLEX CHARGED - SECOND SUBSTRATE CHARGED
								if((nrgBooleanFunction[0] == TRUENRG) || (nrgBooleanFunction[2] == TRUENRG) || (nrgBooleanFunction[4] == TRUENRG))
								{
									// Compute the total amount of complexes and second substrates
									temp_cpxAmount = speciesIter->getChargeMols();
									temp_substrateAmount = allSpecies.at(temp_mol_II).getChargeMols();
									performSingleGilleSpieIntroduction(temp_cpxAmount, temp_substrateAmount, complexID, temp_mol_II, temp_catalysisID, ENDO_CONDENSATION,
																	   temp_mol_I, temp_mol_II, temp_mol_III, temp_mol_IV, BOTHLOAD, temp_reactionID,false);
								}

								// COMPLEX NOT CHARGED - SECOND SUBSTRATE NOT CHARGED (ONLY IF CONDENSATION ESOERGONIC)
								if(nrgBooleanFunction[7] == TRUENRG)
								{
									// Compute the total amount of complexes and second substrates
									temp_cpxAmount = speciesIter->getNOTchargeMols();
									temp_substrateAmount = allSpecies.at(temp_mol_II).getNOTchargeMols();
									performSingleGilleSpieIntroduction(temp_cpxAmount, temp_substrateAmount, complexID, temp_mol_II, temp_catalysisID, CONDENSATION,
																	   temp_mol_I, temp_mol_II, temp_mol_III, temp_mol_IV, NOTHINGLOAD, temp_reactionID,false);
								}
							}else{ // NO ENERGY
								// Compute the total amount of complexes and second substrates
								temp_cpxAmount = speciesIter->getAmount();
								temp_substrateAmount_TOT = allSpecies.at(temp_mol_II).getAmount();
								performSingleGilleSpieIntroduction(temp_cpxAmount, temp_substrateAmount_TOT, complexID, temp_mol_II, temp_catalysisID, CONDENSATION,
																   temp_mol_I, temp_mol_II, temp_mol_III, temp_mol_IV, NOTHINGLOAD, temp_reactionID, false);
							} // end if(nrgBoolFlag == ENERGYBASED)
						}catch(exception&e){
						     cout << "Source Code Line: " << __LINE__ << endl;
						     cout << "CpxAmout -> " << temp_cpxAmount << " - Sub Amount -> " << temp_substrateAmount_TOT << endl;
						     cerr << "exceptioncaught:" << e.what() << endl;
						     ExitWithError("performOPTGillespieComputation::2nd step condensation","exceptionerrorthrown");
						}
					} // for(acs_int listCondSecStep = 0; listCondSecStep < speciesIter->getSecSubListSize(); listCondSecStep++)
				} // end if(allSpecies.at(mid).getSecSubListSize() > 0)
			} //end if total amount greater than 0
		} // end for each species

		// IF THERE ARE CATALYSIS ************************************************************************************

		if((acs_longInt)allCatalysis.size() > 0)
		{
			for(vector<catalysis>::iterator catalysisIter = allCatalysis.begin(); catalysisIter != allCatalysis.end(); catalysisIter++)
			{
				tempIdCatalysis = catalysisIter->getCatId();
				// Retrive the energetic boolean function of the reaction
				try{
					nrgBooleanFunction = dec2bin(allReactions.at(catalysisIter->getReactionID()).getEnergyType());
				}catch(exception&e){
					cout << "nrgBooleanFunction = dec2bin(allReactions.at(allCatalysis.at(idCat).getReactionID()).getEnergyType());" << endl;
					cout << "Vectorsize "<<allReactions.size() << " - position " << catalysisIter->getReactionID() << endl;
					cout << "Line -> " << __LINE__ << endl;
					cerr << "exceptioncaught:" << e.what() << endl;
					ExitWithError("performGillespieComputation","exceptionerrorthrown");
			   }

		//\-- CLEAVAGE .-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

				//IF CLEAVAGE REACTION (if reverse reactions are enabled reaction is computed as well)
			   // TRY HAS BEEN ALREADY USED
				if((allReactions.at(catalysisIter->getReactionID()).getType() == CLEAVAGE) || (reverseReactions == true))
				{
					temp_mol_I = allReactions.at(catalysisIter->getReactionID()).getSpecies_I();     // SUBSTRATE
					temp_mol_II = allReactions.at(catalysisIter->getReactionID()).getSpecies_II();   // PRODUCT 1
					temp_mol_III = allReactions.at(catalysisIter->getReactionID()).getSpecies_III(); // PRODUCT 2
					temp_mol_IV = catalysisIter->getCat();	// CATALYST
					temp_reactionID = catalysisIter->getReactionID(); // REACTION ID

					// Compute the overall number of molecules for the species involved
					temp_catAmount_TOT = allSpecies.at(temp_mol_IV).getAmount();
					temp_catAmount_charged = allSpecies.at(temp_mol_IV).getChargeMols();
					temp_catAmount_NotCharged = allSpecies.at(temp_mol_IV).getNOTchargeMols();

				   try{
						temp_substrateAmount_TOT = allSpecies.at(temp_mol_I).getAmount();
						temp_substrateAmount_charged = allSpecies.at(temp_mol_I).getChargeMols();
						temp_substrateAmount_NotCharged = allSpecies.at(temp_mol_I).getNOTchargeMols();
					}catch(exception&e){
						cout << "temp_substrateAmount_TOT = allSpecies.at(temp_mol_I).getAmount();" << endl;
						cout << "Vectorsize "<< allSpecies.size() << " - position " << temp_mol_I << endl;
						cerr << "exceptioncaught:" << e.what() << endl;
						ExitWithError("performGillespieComputation","exceptionerrorthrown");
				   }

					// Gillespie record creation according to the energy configuration
					try{
						if(nrgBoolFlag == ENERGYBASED)
						{

							// CATALYST NOT LOADED, SUBSTRATE LOADED (- -)
							if(nrgBooleanFunction[11] == TRUENRG)
							if(checkAvailability(temp_mol_IV, temp_mol_I, temp_catAmount_NotCharged, temp_substrateAmount_NotCharged))
								performSingleGilleSpieIntroduction(temp_catAmount_NotCharged, temp_substrateAmount_NotCharged, temp_mol_IV, temp_mol_I, tempIdCatalysis, CLEAVAGE,
									   temp_mol_I, temp_mol_II, temp_mol_III, temp_mol_IV, NOTHINGLOAD, temp_reactionID, true);

							// CATALYST NOT LOADED, SUBSTRATE LOADED (- +)
							if(nrgBooleanFunction[10] == TRUENRG)
									performSingleGilleSpieIntroduction(temp_catAmount_NotCharged, temp_substrateAmount_charged, temp_mol_IV, temp_mol_I, tempIdCatalysis, ENDO_CLEAVAGE,
																		temp_mol_I, temp_mol_II, temp_mol_III, temp_mol_IV, SUBSTRATELOAD, temp_reactionID, false);
							// CATALYST NOT LOADED, SUBSTRATE LOADED (+ -)
							if(nrgBooleanFunction[9] == TRUENRG)
									performSingleGilleSpieIntroduction(temp_catAmount_charged, temp_substrateAmount_NotCharged, temp_mol_IV, temp_mol_I, tempIdCatalysis, ENDO_CLEAVAGE,
																	   temp_mol_I, temp_mol_II, temp_mol_III, temp_mol_IV, CATALYSTLOAD, temp_reactionID, false);

							// CATALYST NOT LOADED, SUBSTRATE LOADED (+ +)
							if(nrgBooleanFunction[8] == TRUENRG)
								if(checkAvailability(temp_mol_IV, temp_mol_I, temp_catAmount_charged, temp_substrateAmount_charged))
									performSingleGilleSpieIntroduction(temp_catAmount_charged, temp_substrateAmount_charged, temp_mol_IV, temp_mol_I, tempIdCatalysis, ENDO_CLEAVAGE,
									 temp_mol_I, temp_mol_II, temp_mol_III, temp_mol_IV, BOTHLOAD, temp_reactionID, true);
						}else{
								// CATALYST NOT LOADED, SUBSTRATE NOT LOADED
							if(checkAvailability(temp_mol_IV, temp_mol_I, temp_catAmount_TOT, temp_substrateAmount_TOT))
								performSingleGilleSpieIntroduction(temp_catAmount_TOT, temp_substrateAmount_TOT, temp_mol_IV, temp_mol_I, tempIdCatalysis, CLEAVAGE,
																   temp_mol_I, temp_mol_II, temp_mol_III, temp_mol_IV, NOTHINGLOAD, temp_reactionID, true);
						} // end if(nrgBoolFlag == ENERGYBASED)
					}catch(exception&e){
					     cout << "Source Code Line: " << __LINE__ << endl;
					     cerr << "exceptioncaught:" << e.what() << endl;
					     ExitWithError("performOPTGillespieComputation::Cleavage","exceptionerrorthrown");
					}
				} // end getType() == CLEAVAGE

		//\-- OVERALL CONDENSATION .-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-

				// IF THIS IS A CONDENSATION REACTION (if reverse reactions are enabled reaction is comuputed as well)
				if((allReactions.at(catalysisIter->getReactionID()).getType() == CONDENSATION) || (reverseReactions == true))
				{

					//\-- COMPLEX CREATION	.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

					temp_mol_I = catalysisIter->getCat(); // Catalyst
					// According to the complex creation target, species to create the complex is selected.
					if(catalysisIter->getCpxTarget() == 1)
					{
						temp_mol_II = allReactions.at(catalysisIter->getReactionID()).getSpecies_II(); // First substrate
						temp_mol_IV = allReactions.at(catalysisIter->getReactionID()).getSpecies_III(); // Second substrate
					}else{
						temp_mol_II = allReactions.at(catalysisIter->getReactionID()).getSpecies_III(); // First substrate
						temp_mol_IV = allReactions.at(catalysisIter->getReactionID()).getSpecies_II(); // Second substrate
					}

					//temp_mol_III = 0; // This variable contains the number of energized molecules (0 ++ 2 +- 4 -+ 6 --)
					temp_mol_III = catalysisIter->getCatId();; // This variable contain the catalysis ID
					temp_reactionID = catalysisIter->getReactionID();

					//\-- FIRST AND SECOND SUBSTRATE  -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-

					// Compute the overall number of molecules for the species involved
					temp_catAmount_TOT = allSpecies.at(temp_mol_I).getAmount();
					temp_catAmount_charged = allSpecies.at(temp_mol_I).getChargeMols();
					temp_catAmount_NotCharged = allSpecies.at(temp_mol_I).getNOTchargeMols();

					temp_substrateAmount_TOT = allSpecies.at(temp_mol_II).getAmount();
					temp_substrateAmount_charged = allSpecies.at(temp_mol_II).getChargeMols();
					temp_substrateAmount_NotCharged = allSpecies.at(temp_mol_II).getNOTchargeMols();

					try{

						if(nrgBoolFlag == ENERGYBASED)
						{

							// CATALYST LOADED, SUBSTRATE LOADED (+ + ...)
							if((nrgBooleanFunction[0] == TRUENRG) || (nrgBooleanFunction[1] == TRUENRG))
							if(checkAvailability(temp_mol_I, temp_mol_II, temp_catAmount_charged, temp_substrateAmount_charged))
								performSingleGilleSpieIntroduction(temp_catAmount_charged, temp_substrateAmount_charged, temp_mol_I, temp_mol_II, tempIdCatalysis, ENDO_COMPLEXFORMATION,
										temp_mol_I, temp_mol_II, temp_mol_III, temp_mol_IV, BOTHLOAD, temp_reactionID, true);

							// CATALYST LOADED, SUBSTRATE NOT LOADED (+ - ...)
							if((nrgBooleanFunction[2] == TRUENRG) || (nrgBooleanFunction[3] == TRUENRG))
									performSingleGilleSpieIntroduction(temp_catAmount_charged, temp_substrateAmount_NotCharged, temp_mol_I, temp_mol_II, tempIdCatalysis, ENDO_COMPLEXFORMATION,
																		temp_mol_I, temp_mol_II, temp_mol_III, temp_mol_IV, CATALYSTLOAD, temp_reactionID, false);

							// CATALYST NOT LOADED, SUBSTRATE LOADED (- + ...)
							if((nrgBooleanFunction[4] == TRUENRG) || (nrgBooleanFunction[5] == TRUENRG))
									performSingleGilleSpieIntroduction(temp_catAmount_NotCharged, temp_substrateAmount_charged, temp_mol_I, temp_mol_II, tempIdCatalysis, ENDO_COMPLEXFORMATION,
																		temp_mol_I, temp_mol_II, temp_mol_III, temp_mol_IV, SUBSTRATELOAD, temp_reactionID, false);

							// CATALYST NOT LOADED, SUBSTRATE NOT LOADED (- - ...)
							if((nrgBooleanFunction[6] == TRUENRG) || (nrgBooleanFunction[7] == TRUENRG))
								if(checkAvailability(temp_mol_I, temp_mol_II, temp_catAmount_NotCharged, temp_substrateAmount_NotCharged))
									performSingleGilleSpieIntroduction(temp_catAmount_NotCharged, temp_substrateAmount_NotCharged, temp_mol_I, temp_mol_II, tempIdCatalysis, COMPLEXFORMATION,
										temp_mol_I, temp_mol_II, temp_mol_III, temp_mol_IV, NOTHINGLOAD, temp_reactionID,true);
						}else{ // NO ENERGY
						   if(checkAvailability(temp_mol_I, temp_mol_II, temp_catAmount_TOT, temp_substrateAmount_TOT))
								performSingleGilleSpieIntroduction(temp_catAmount_TOT, temp_substrateAmount_TOT, temp_mol_I, temp_mol_II, tempIdCatalysis, COMPLEXFORMATION,
									   temp_mol_I, temp_mol_II, temp_mol_III, temp_mol_IV, NOTHINGLOAD, temp_reactionID,true);
						} // end if(nrgBoolFlag == ENERGYBASED)
					}catch(exception&e){
					     cout << "Source Code Line: " << __LINE__ << endl;
					     cerr << "exceptioncaught:" << e.what() << endl;
					     ExitWithError("performOPTGillespieComputation::1nd step condensation","exceptionerrorthrown");
					}
				} // end if(allReactions.at(allCatalysis.at(idCat).getReactionID()).getType() == CONDENSATION)
			} // end for all catalysis
		} // end if if((acs_longInt)allCatalysis.size() > 0)

		// -----------------------------------------------------------------------
		// SPONTANEOUS REACTIONS SECTION
		// IF SPONTANEOUS REACTIONS ARE TURNED ON, ALL REACTIONS MUST BE EVALUATED
		// -----------------------------------------------------------------------

		if(spontRct > 0)
		{
			for(vector<reactions>::iterator reactionsIter = allReactions.begin(); reactionsIter != allReactions.end(); reactionsIter++)
			{
				// If the spontanoues constant of the reaction is greater than 0
				if(reactionsIter->getKspont() > 0)
				{
					temp_mol_I = reactionsIter->getSpecies_I(); // Mol_I
					temp_mol_II = reactionsIter->getSpecies_II(); // Mol_II
					temp_mol_III = reactionsIter->getSpecies_III(); // Mol_III
					temp_mol_IV = 0; // To NOT confuse with species ID
					temp_catalysisID = 0; // To NOT confuse with catalysis 0, no catalysis are involved in this process.
					temp_reactionID = reactionsIter->getID(); // reaction ID
					temp_rctType = reactionsIter->getType();; // reaction type
					if((temp_rctType == CLEAVAGE) || (temp_rctType == ENDO_CLEAVAGE))
					{
						acs_double tempScore = allSpecies.at(temp_mol_I).getAmount() * reactionsIter->getKspont(); // uni-molecular reaction score in case of cleavage
						if(debugLevel >= MEDIUM_DEBUG) cout << "\t\t\t|- SPONTANEOUS DISSOCIATION: " << tempScore << endl;
						if(tempScore > 0)
						{
							allGillespieScores.push_back(gillespie((acs_longInt)allGillespieScores.size(), SPONTANEOUS_CLEAVAGE, tempScore,
														 temp_mol_I, temp_mol_II, temp_mol_III,temp_mol_IV, temp_reactionID, temp_catalysisID));
							gillespieTotalScore += tempScore;
							gillespieCumulativeStepScoreList.push_back(gillespieTotalScore);
							// UPDATE SPECIES GILLESPIE ENGAGEMENT
							allSpecies.at(temp_mol_I).insertGillID(allGillespieScores.back().getID());
							// In the case of cleavage molII and molIII are products, if they are not evaluated yet, hence probability of new species increases
							if((allSpecies.at(temp_mol_II).getEvaluated() == 0) || (allSpecies.at(temp_mol_III).getEvaluated() == 0)) gillespieNewSpeciesScore += tempScore;
						}

					}
					try{
						if((temp_rctType == CONDENSATION) || (temp_rctType == ENDO_CONDENSATION))
						{ // Spontaneous condensation is a bi-molecolar reaction
							if(debugLevel >= MEDIUM_DEBUG) cout << "\t\t\t|- SPONTANEOUS CONDENSATION" << endl;
							if(checkAvailability(temp_mol_II, temp_mol_III, allSpecies.at(temp_mol_II).getAmount(), allSpecies.at(temp_mol_III).getAmount()))
								performSingleGilleSpieIntroduction(allSpecies.at(temp_mol_II).getAmount(), allSpecies.at(temp_mol_III).getAmount(),
																   temp_mol_II, temp_mol_III, tempIdCatalysis, SPONTANEOUS_CONDENSATION,
																   temp_mol_I, temp_mol_II, temp_mol_III, temp_mol_IV, NOTHINGLOAD, temp_reactionID, true);
						}
					}catch(exception&e){
					     cout << "Source Code Line: " << __LINE__ << endl;
					     cerr << "exceptioncaught:" << e.what() << endl;
					     ExitWithError("performOPTGillespieComputation::spontaneous","exceptionerrorthrown");
					}
				}
			}
			if(debugLevel == SMALL_DEBUG) cout << "step " << tmpActSTEP << endl;
			if(debugLevel == SMALL_DEBUG) printGillespieStructure();
		} // end IF SPONTANEOUS REACTIONS ARE TURNED ON, ALL REACTIONS MUST BE EVALUATED
		if(debugLevel == SMALL_DEBUG)
		{
			cout << "\t\t|- GILLESPIE STRUCTURE CREATION ENDED. Number of possible reactions:" << allGillespieScores.size() << endl;
		}
    } // end if(tmpActSTEP == 1)


    // Store gillespie computational time and start performReaction Times
    gillespiePartialTimes.push_back(((float)clock() - gillespiePartialTimer) / CLOCKS_PER_SEC);
    performReactionPartialTimer = clock();

	bool goReaction = true;
	acs_double minimalTimeForOneMols = 1 / (influx_rate * AVO);
	if((acs_longInt)allGillespieScores.size() > 0)
	{
            if(debugLevel == SMALL_DEBUG) printGillespieStructure();

            // SELECT REACTION WITHIN THE GILLESPIE STRUCTURE
            if((acs_longInt)allGillespieScores.size() == 1)
            {
                reaction_u = 0;
            }else{
            	try{
					reaction_u = returnSelectionIdFromAWeightProbVector(gillespieCumulativeStepScoreList, gillespieCumulativeStepScoreList.back(), tmpRndDoubleGen, __LINE__);
					if(debugLevel == GILLESPIESTUFF){
						printGillespieStructure();
						printInitialCondition();
						cin.ignore().get();
						//cout <<  allGillespieScores.at(reaction_u) << endl;
					}
            	}catch(exception&e){
				     cout << "Source Code Line: " << __LINE__ << endl;
				     cerr << "exceptioncaught:" << e.what() << endl;
				     ExitWithError("performOPTGillespieComputation::Gillespie selection","exceptionerrorthrown");
				}
            }

            // CREATE RANDOM NUMBER TO COMPUTE THE TIME

            tmpDeltaT = ((1 / gillespieTotalScore) * log(1 / tmpRndDoubleGen()));

            // If deltaT is lower than 10 seconds, it is fixed to 10 second in order to continue the simulation
            if((tmpDeltaT) > MINIMALRCTTIMEMULTI*minimalTimeForOneMols)
            {
				tmpDeltaT = MINIMALRCTTIMEMULTI*minimalTimeForOneMols;
				goReaction = false;
            }

            if(debugLevel == SMALL_DEBUG)
                cout << "gillespieTotalScore: " << gillespieTotalScore << " G size: " << allGillespieScores.size() << " delta T " << tmpDeltaT << endl;

            if(goReaction)
            {
                if(reactionsTime.size() > 0)
                {
                    tempTime = reactionsTime.at(reactionsTime.size() - 1) + tmpDeltaT;
                    timeSinceTheLastInFlux += tmpDeltaT;
                }else{
                    tempTime = 0.0;
                }
                reactionsTime.push_back(tempTime);
                setActualTime(tempTime);
                gillespieReactionsSelected.push_back(reaction_u);
                allTimes.push_back(((float)clock() - tmpTimeElapsed) / CLOCKS_PER_SEC);

                // =^=^=^=^=^=^=^=^=^=^=^=^=^=^=^=^=^=^=^=^=^=^=^=^
                // PERFORM REACTION SELECTED BEFORE
                // ^=^=^=^=^=^=^=^=^=^=^=^=^=^=^=^=^=^=^=^=^=^=^=^=
                // Compute Gillespie mean
                gillespieMean = gillespieTotalScore / (acs_longInt)allGillespieScores.size();
                if(gillespieTotalScore > 0){ratioBetweenNewGillTotGill = gillespieNewSpeciesScore / gillespieTotalScore;}else{ratioBetweenNewGillTotGill=0;}
                if(gillespieTotalScore > 0){ratioBetweenReverseAndTotalScore = reverseReactionsGillScore / gillespieTotalScore;}else{ratioBetweenReverseAndTotalScore=0;}
                if(!devStd()) // compute Gillespie score vector standard deviation
                     ExitWithError("devStd", "Problems during Gillespie score standard deviation computation");
                if(!entropy()) // compute Gillespie score vector entropy
                     ExitWithError("entropy", "Problems during Gillespie score entropy computation");

                // PERFORM REACTION
                try{
                	if(!performReaction(reaction_u, tmpRndDoubleGen, tmpActGEN, tmpActSIM, tmpActSTEP, tmpStoringPath))
                                ExitWithError("performReaction", "Problems during the reaction computation");
                }catch(exception&e){
				     cout << "Source Code Line: " << __LINE__ << endl;
				     cerr << "exceptioncaught:" << e.what() << endl;
				     ExitWithError("performOPTGillespieComputation::Perform Reaction","exceptionerrorthrown");
				}

                // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                // CHANGE VOLUME IF PROTOCELL WITH VARYING VOLUME
                // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                if(volumeGrowth) changeVolume(tmpDeltaT);

            }else{
				gillespieMean = 0;
				gillespieSD = 0;
				gillespieEntropy = 0;
				if(reactionsTime.size() > 0)
				{
					tempTime = reactionsTime.at(reactionsTime.size() - 1) + MINIMALRCTTIMEMULTI*minimalTimeForOneMols;
					timeSinceTheLastInFlux += MINIMALRCTTIMEMULTI*minimalTimeForOneMols;
					tmpDeltaT = MINIMALRCTTIMEMULTI*minimalTimeForOneMols;
				}else{
					tempTime = 0.0;
				}
				reactionsTime.push_back(tempTime);
				setActualTime(tempTime);
				gillespieReactionsSelected.push_back(0);
				allTimes.push_back(((float)clock() - tmpTimeElapsed) / CLOCKS_PER_SEC);

				if(debugLevel >= RUNNING_VERSION)
								cout << "\t\t\t|- NO REACTIONS AT THIS STEP T:" << tempTime << " G: " << allGillespieScores.size() << endl;
            }

        }else{ // If there are not possible reactions
            gillespieMean = 0;
            gillespieSD = 0;
            gillespieEntropy = 0;
            ratioBetweenNewGillTotGill = 0;
            if(reactionsTime.size() > 0)
            {
                    tempTime = reactionsTime.at(reactionsTime.size() - 1) + MINIMALRCTTIMEMULTI*minimalTimeForOneMols;
                    timeSinceTheLastInFlux += MINIMALRCTTIMEMULTI*minimalTimeForOneMols;
                    tmpDeltaT = MINIMALRCTTIMEMULTI*minimalTimeForOneMols;
            }else{
                    tempTime = MINIMALRCTTIMEMULTI*minimalTimeForOneMols;
                    timeSinceTheLastInFlux = MINIMALRCTTIMEMULTI*minimalTimeForOneMols;
                    tmpDeltaT = MINIMALRCTTIMEMULTI*minimalTimeForOneMols;
            }
            reactionsTime.push_back(tempTime);
            setActualTime(tempTime);
            gillespieReactionsSelected.push_back(0);
            allTimes.push_back(((float)clock() - tmpTimeElapsed) / CLOCKS_PER_SEC);

            if(debugLevel >= RUNNING_VERSION)
                    cout << "\t\t\t|- NO REACTIONS AT THIS STEP" << endl;
        } // end if((acs_longInt)allGillespieScores.size() > 0)

        // Store perform reaction time and start remaining processes timer
        performReactionPartialTimes.push_back(((float)clock() - performReactionPartialTimer) / CLOCKS_PER_SEC);
        remainingProcessesPartialTimer = clock();

        // If the system is open influx and efflux processes are performed
        if(influx_rate > 0)
        {
        	//     acs_double minimalTimeForOneMols = 1 / (influx_rate * AVO);
			if(debugLevel >= SMALL_DEBUG)
			{
				cout << "\t\t\t|- REFILL --------" << endl;
				cout << "\t\t\t\t|- Minimal Time for one molecule: " << minimalTimeForOneMols <<
												" - Time since the last influx: " << timeSinceTheLastInFlux <<  endl;
			}
			// If the time interval betweem two successive influx is enough to introduce at least one new mol...
			if(timeSinceTheLastInFlux > minimalTimeForOneMols)
			{
				if(debugLevel >= SMALL_DEBUG)
				{
					cout << "\t\t\t\t|- Time: " << reactionsTime.at(reactionsTime.size() - 1)
						 << " - Time needed for 1 molecule incoming: " << minimalTimeForOneMols
						 << " - Time Since The Last InFlux: " << timeSinceTheLastInFlux << endl;
				}
				// PERFORM REFILL !!!
				try{
					if(performRefill(timeSinceTheLastInFlux, minimalTimeForOneMols, tmpRndDoubleGen))
									timeSinceTheLastInFlux = 0;
				}catch(exception&e){
				     cout << "Source Code Line: " << __LINE__ << endl;
				     cerr << "exceptioncaught:" << e.what() << endl;
				     ExitWithError("performOPTGillespieComputation::Perform Refill","exceptionerrorthrown");
				}
			}
			// PERFORM EFFLUX PROCESS
			try{
				performMoleculesEfflux(tmpDeltaT, tmpRndDoubleGen);
			}catch(exception&e){
			     cout << "Source Code Line: " << __LINE__ << endl;
			     cerr << "exceptioncaught:" << e.what() << endl;
			     ExitWithError("performOPTGillespieComputation::Perform molecules efflux","exceptionerrorthrown");
			}
         }

    // Perform molecule charging
    if(energy == ENERGYBASED)
        performDETMoleculesCharging(tmpDeltaT, tmpRndDoubleGen);

    try{
		// Perform deterministic complex dissociation
		if(!performDETComplexDissociation(tmpDeltaT, tmpRndDoubleGen))
			ExitWithError("performComplexDissociation", "Problems during the COMPLEX DETERMINISTIC DISSOCIATION process in Gillespie");
    }catch(exception&e){
	     cout << "Source Code Line: " << __LINE__ << endl;
	     cerr << "exceptioncaught:" << e.what() << endl;
	     ExitWithError("performOPTGillespieComputation::complex dissociation","exceptionerrorthrown");
	}

    // COPY and CLEAN OF GILLESPIE STRUCTURES
    COPYOFallGillespieScores = allGillespieScores;
    allGillespieScores.clear();
    gillespieCumulativeStepScoreList.clear();

    // Store remaining processes time
    remainingProcessesPartialTimes.push_back(((float)clock() - remainingProcessesPartialTimer) / CLOCKS_PER_SEC);

    if(debugLevel == FINDERRORDURINGRUNTIME) cout << "environment::performOPTGillespieComputation end" << endl;

    return flagControl;

}


/**
 Compute a single gillespie score according to the amount and peoprieties of the species involved
 @version 1.0
 @date 20110214
 */
acs_double environment::computeSinglGilScore(acs_longInt tmpAmountI, acs_double tmpDifI, acs_int tmpSolI,
                                             acs_longInt tmpAmountII, acs_double tmpDifII, acs_int tmpSolII, acs_double tmpK, bool tmpSameMol)
{
	if(debugLevel == FINDERRORDURINGRUNTIME) cout << "\t\tenvironment::computeSinglGilScore start" << endl;
	
	acs_double tempScore;
	
	if(tmpSameMol) //IF THE SPECIES ARE THE SAME
	{
		tempScore = ((tmpAmountI * tmpDifI * tmpSolI * (tmpAmountII - 1) * tmpDifII * tmpSolII) / 2) * (tmpK / (AVO * volume));
	}else{
		tempScore = tmpAmountI * tmpDifI * tmpSolI * tmpAmountII * tmpDifII * tmpSolII * tmpK / (AVO * volume);
	}
	
	// cia
	
	if(debugLevel == FINDERRORDURINGRUNTIME) cout << "\t\tenvironment::computeSinglGilScore end" << endl;
	
	return tempScore;
}

/**
 Compute and introduct a single Gillespie entry within the Gillespie Structure
 @version 1.0
 @date 20110222
 */
void environment::performSingleGilleSpieIntroduction(acs_longInt tmpAmountI, acs_longInt tmpAmountII, acs_longInt tmpIDI, acs_longInt tmpIDII, acs_longInt tmpIDCatalysis,
													 acs_int tmp__rctType, acs_longInt tmpMol_I, acs_longInt tmpMol_II, acs_longInt tmpMol_III, acs_longInt tmpMol_IV,
													 acs_int tmp_NRGDirection, acs_longInt tmpRctID, bool tmpSameSpeciesControl)
{
    if(debugLevel == FINDERRORDURINGRUNTIME) cout << "\tenvironment::performSingleGilleSpieIntroduction start" << endl;
    bool temp_sameSpecies = false;
    acs_double temp_score; // Gillespie score
    acs_double temp_k_reaction; // k reaction

    if((tmpAmountI > 0) && (tmpAmountII > 0))
    {
        try{
            switch(tmp__rctType)
            {
                case CONDENSATION:
					temp_k_reaction = allCatalysis.at(tmpIDCatalysis).getKass();
					break;
                case ENDO_CONDENSATION:
					temp_k_reaction = allCatalysis.at(tmpIDCatalysis).getKass();
					break;
                case CLEAVAGE:
					temp_k_reaction = allCatalysis.at(tmpIDCatalysis).getKdiss();
					break;
                case ENDO_CLEAVAGE:
					temp_k_reaction = allCatalysis.at(tmpIDCatalysis).getKdiss();
					break;
                case COMPLEXFORMATION:
					temp_k_reaction = allCatalysis.at(tmpIDCatalysis).getK_cpx();
					break;
                case ENDO_COMPLEXFORMATION:
                    temp_k_reaction = allCatalysis.at(tmpIDCatalysis).getK_cpx();
					break;
                case SPONTANEOUS_CONDENSATION:
                	temp_k_reaction = allReactions.at(tmpRctID).getKspont();
                	break;
                case SPONTANEOUS_CLEAVAGE:
                	temp_k_reaction = allReactions.at(tmpRctID).getKspont();
            }
            if(tmpSameSpeciesControl)
            {
                if(tmpIDI == tmpIDII)
                        temp_sameSpecies = true;
            }
        }catch(exception&e){
            cout << "switch(tmp__rctType)" << endl;
            cout << "Vectorsize "<< allCatalysis.size() << " - position " << tmpIDCatalysis << endl;
            cerr << "exceptioncaught:" << e.what() << endl;
            ExitWithError("performSimgleGilleSpieIntroduction","exceptionerrorthrown");
       }

        try{
            temp_score = computeSinglGilScore(tmpAmountI, allSpecies.at(tmpIDI).getDiffusionEnh(), allSpecies.at(tmpIDI).getSolubility(),
                                              tmpAmountII, allSpecies.at(tmpIDII).getDiffusionEnh(), allSpecies.at(tmpIDII).getSolubility(),
                                              temp_k_reaction, temp_sameSpecies);

        }catch(exception&e){
            cout << "temp_score = computeSinglGilScore" << endl;
            cout << "Vectorsize "<< allSpecies.size() << " - position " << tmpIDI << endl;
            cerr << "exceptioncaught:" << e.what() << endl;
            ExitWithError("performSimgleGilleSpieIntroduction","exceptionerrorthrown");
       }

        if(temp_score > 0) // If score is grater than 0
        {
            if((debugLevel >= SMALL_DEBUG))
                    cout << "\t\t\t|- Amount 1: " << tmpIDI << " " << tmpAmountI <<
                    " rct type: " << tmp__rctType <<
                    " Amount 2: " << tmpIDII << " "  << tmpAmountII <<
                    " Reaction: " << tmpRctID <<
                    " tmp_NRGDirection: " << tmp_NRGDirection <<
                    " tmpIDCatalysis: " << tmpIDCatalysis <<
                    " cpxTarget: " << allCatalysis.at(tmpIDCatalysis).getCpxTarget() <<
                    " cat: " << allCatalysis.at(tmpIDCatalysis).getCat() <<
                    " mol1: " << allReactions.at(allCatalysis.at(tmpIDCatalysis).getReactionID()).getSpecies_I() <<
                    " mol2: " << allReactions.at(allCatalysis.at(tmpIDCatalysis).getReactionID()).getSpecies_II() <<
                    " mol3: " << allReactions.at(allCatalysis.at(tmpIDCatalysis).getReactionID()).getSpecies_III() << endl;

            // Check whether the gillespie new entry is already present, if not new entry will be added
			try{
				bool gillAlreadyPresent = false;
				if(allGillespieScores.size() > 0)
				{
					for(vector<gillespie>::iterator tmpAllGillIter = allGillespieScores.begin(); tmpAllGillIter != allGillespieScores.end(); tmpAllGillIter++)
					{
						if((tmpAllGillIter->getMolI() == tmpMol_I) &&
						   (tmpAllGillIter->getMolII() == tmpMol_II) &&
						   (tmpAllGillIter->getMolIII() == tmpMol_III) &&
						   (tmpAllGillIter->getMolIV() == tmpMol_IV) &&
						   (tmpAllGillIter->getIdReactionType() == tmp__rctType))
						{
							gillAlreadyPresent = true;
							break;
						}
					}

				}

				if(!gillAlreadyPresent)
				{
					try{
						//cout << tmpRctID << " " << tmpIDCatalysis << " " << tmpMol_I << " " << tmpMol_II<< endl;
						allGillespieScores.push_back(gillespie((acs_longInt)allGillespieScores.size(),
													tmp__rctType, temp_score, tmpMol_I, tmpMol_II,
													tmpMol_III, tmpMol_IV, tmp_NRGDirection, tmpRctID,
													tmpIDCatalysis));
						gillespieTotalScore += temp_score;
						gillespieCumulativeStepScoreList.push_back(gillespieTotalScore);
						// Update single species gillespie engagements list
/*						try{ // Until the new Gillespie computation enhancement these lines are commented.
							allSpecies.at(tmpIDI).insertGillID(allGillespieScores.back().getID());
							allSpecies.at(tmpIDII).insertGillID(allGillespieScores.back().getID());
						}catch(exception&e){
							cout << "!!! ERROR in line " << __LINE__ << endl;
							cout << "|- Gill size -> " << allGillespieScores.size() << " || gillCumuStepScoreList size -> " << gillespieCumulativeStepScoreList.size() << endl
								 << " || tmpIDI: " << tmpIDI << " || tmpIDII: " << tmpIDII << endl
								 << " || M1: " << tmpMol_I << " || M2: " << tmpMol_II << " || M3: " << tmpMol_III << " || M4: " << tmpMol_IV
								 << " || type: " << tmp__rctType << " || RCTID: "  << tmpRctID << " || Catalysis: " << tmpIDCatalysis << endl;
							cerr << "exceptioncaught:" << e.what() << endl;
							ExitWithError("performSimgleGilleSpieIntroduction","exceptionerrorthrown");
						}*/
							// If the theoretical product is not evaluated gillespieNewSpeciesScore is incremented
						if((tmp__rctType == CONDENSATION) || (tmp__rctType == ENDO_CONDENSATION) || (tmp__rctType == SPONTANEOUS_CONDENSATION))
						{
							// IN the case of condensation molIII is the product,
							// in the reaction structure molIII is a substrate but this subroutine is the product,
							// mol_I is the complex and mol_II is the second substrate
							if(allSpecies.at(tmpMol_III).getEvaluated() == 0)
								gillespieNewSpeciesScore += temp_score;
						}else if((tmp__rctType == CLEAVAGE) || (tmp__rctType == ENDO_CLEAVAGE) || (tmp__rctType == SPONTANEOUS_CLEAVAGE))
						{
							// In the case of cleavage molII and molIII are products
							if((allSpecies.at(tmpMol_II).getEvaluated() == 0) || (allSpecies.at(tmpMol_III).getEvaluated() == 0))
								gillespieNewSpeciesScore += temp_score;
						}

						if( ((tmp__rctType == CLEAVAGE || tmp__rctType == ENDO_CLEAVAGE) && (allReactions.at(tmpRctID).getType() == CONDENSATION))
							 ||
							((tmp__rctType == CONDENSATION || tmp__rctType == ENDO_CONDENSATION) && (allReactions.at(tmpRctID).getType() == CLEAVAGE)) )
						{
							if(reverseReactions == false)
							{
								printGillespieStructure();
								ExitWithError("performSimgleGilleSpieIntroduction","Reverse reaction is detected although it shouldn't be possible...");
							}
							reverseReactionsGillScore += temp_score;

						}
					}catch(exception&e){
		                cout << "!!! ERROR in line " << __LINE__ << endl;
		                cout << "|- Gill size -> " << allGillespieScores.size() << " || gillCumuStepScoreList size -> " << gillespieCumulativeStepScoreList.size()
		                	 << " || M1: " << tmpMol_I << " || M2: " << tmpMol_II << " || M3: " << tmpMol_III << " || M4: " << tmpMol_IV
		                	 << " || type: " << tmp__rctType << " || RCTID: "  << tmpRctID << " || Catalysis: " << tmpIDCatalysis << endl;
		                cerr << "exceptioncaught:" << e.what() << endl;
		                ExitWithError("performSimgleGilleSpieIntroduction","exceptionerrorthrown");
		           }
				}
			}catch(exception&e){
                cout << "!!! ERROR in line " << __LINE__ << endl;
                cerr << "exceptioncaught:" << e.what() << endl;
                ExitWithError("performSimgleGilleSpieIntroduction","exceptionerrorthrown");
           }
        }
    }
    if(debugLevel == FINDERRORDURINGRUNTIME) cout << "\tenvironment::performSingleGilleSpieIntroduction end" << endl;
}


/**
 This Function returns the total amount of MOLECULES in the  tmpSpeciesvector
 The file is saved in the directory indicated as a second parameter in the run command
 @param vector<species>* tmpSpeciesvector pointer to tmpSpeciesvector
 */
acs_longInt environment::getTotalNumberOfSpecies()
{
	if(debugLevel == FINDERRORDURINGRUNTIME) cout << "\t\tenvironment::getTotalNumberOfMolecules start" << endl;
	
	acs_longInt tmpTotNumOfMols = 0;

        for(vector<species>::iterator tmpAllSpeciesIter = allSpecies.begin(); tmpAllSpeciesIter != allSpecies.end(); tmpAllSpeciesIter++)
        {
            if((tmpAllSpeciesIter->getComplexCutPnt() == 0) &&
                (tmpAllSpeciesIter->getAmount() > 0) &&
                (tmpAllSpeciesIter->getSolubility() == SOLUBLE))
                    tmpTotNumOfMols++;
        }
	
	if(debugLevel == FINDERRORDURINGRUNTIME) cout << "\t\tenvironment::getTotalNumberOfMolecules end" << endl;
	
	return tmpTotNumOfMols;
}

/**
 This Function returns the total amount of MOLECULES in the  tmpSpeciesvector
 The file is saved in the directory indicated as a second parameter in the run command
 @param vector<species>* tmpSpeciesvector pointer to tmpSpeciesvector
 */
acs_longInt environment::getTotalNumberOfPossibleCatalysts()
{
	if(debugLevel == FINDERRORDURINGRUNTIME) cout << "\t\tenvironment::getTotalNumberOfMolecules start" << endl;

	acs_longInt tmpTotNumOfMols = 0;

        for(vector<species>::iterator tmpAllSpeciesIter = allSpecies.begin(); tmpAllSpeciesIter != allSpecies.end(); tmpAllSpeciesIter++)
        {
            if((tmpAllSpeciesIter->getComplexCutPnt() == 0) &&
                (tmpAllSpeciesIter->getSolubility() == SOLUBLE))
                    tmpTotNumOfMols++;
        }

	if(debugLevel == FINDERRORDURINGRUNTIME) cout << "\t\tenvironment::getTotalNumberOfMolecules end" << endl;

	return tmpTotNumOfMols;
}

/**
 This Function returns the total amount of MOLECULES in the  tmpSpeciesvector
 The file is saved in the directory indicated as a second parameter in the run command
 @param vector<species>* tmpSpeciesvector pointer to tmpSpeciesvector
 */
acs_longInt environment::getTotalNumberOfMolecules()
{
    if(debugLevel == FINDERRORDURINGRUNTIME) cout << "\t\tenvironment::getTotalNumberOfMolecules start" << endl;

    acs_longInt tmpTotNumOfMols = 0;
    for(vector<species>::iterator tmpAllSpeciesIter = allSpecies.begin(); tmpAllSpeciesIter != allSpecies.end(); tmpAllSpeciesIter++)
    {
        if((tmpAllSpeciesIter->getComplexCutPnt() == 0) &&
            (tmpAllSpeciesIter->getAmount() > 0) &&
            (tmpAllSpeciesIter->getSolubility() == SOLUBLE))
                tmpTotNumOfMols += tmpAllSpeciesIter->getAmount();
    }

    if(debugLevel == FINDERRORDURINGRUNTIME) cout << "\t\tenvironment::getTotalNumberOfMolecules end" << endl;

    return tmpTotNumOfMols;
}


/**
 This Function returns the total amount of COMPLEX SPECIES in the tmpSpeciesvector
 The file is saved in the directory indicated as a second parameter in the run command
 @param vector<species>* tmpSpeciesvector pointer to tmpSpeciesvector
 */
acs_longInt environment::getTotalNumberOfComplexSpecies()
{
	if(debugLevel == FINDERRORDURINGRUNTIME) cout << "\t\tenvironment::getTotalNumberOfComplexSpecies start" << endl;
	
	acs_longInt tmpTotNumOfMols = 0;
        for(vector<species>::iterator tmpAllSpeciesIter = allSpecies.begin(); tmpAllSpeciesIter != allSpecies.end(); tmpAllSpeciesIter++)
        {
            if((tmpAllSpeciesIter->getComplexCutPnt() > 0) && (tmpAllSpeciesIter->getAmount() > 0))
                    tmpTotNumOfMols++;
        }
	
	if(debugLevel == FINDERRORDURINGRUNTIME) cout << "\t\tenvironment::getTotalNumberOfComplexSpecies end" << endl;
	
	return tmpTotNumOfMols;
}

/**
 This Function returns the total amount of MOLECULES in the  tmpSpeciesvector
 The file is saved in the directory indicated as a second parameter in the run command
 @param vector<species>* tmpSpeciesvector pointer to tmpSpeciesvector
 */
acs_longInt environment::getTotalNumberOfComplexes()
{
    if(debugLevel == FINDERRORDURINGRUNTIME) cout << "\t\tenvironment::getTotalNumberOfComplexes start" << endl;

    acs_longInt tmpTotNumOfMols = 0;
    for(vector<species>::iterator tmpAllSpeciesIter = allSpecies.begin(); tmpAllSpeciesIter != allSpecies.end(); tmpAllSpeciesIter++)
    {
        if((tmpAllSpeciesIter->getComplexCutPnt() > 0) && (tmpAllSpeciesIter->getAmount() > 0))
            tmpTotNumOfMols += tmpAllSpeciesIter->getAmount();
    }

    if(debugLevel == FINDERRORDURINGRUNTIME) cout << "\t\tenvironment::getTotalNumberOfComplexes end" << endl;

    return tmpTotNumOfMols;
}

/**
 This Function returns the total amount of MONOMERS in the  tmpSpeciesvector
 The file is saved in the directory indicated as a second parameter in the run command
 @param vector<species>* tmpSpeciesvector pointer to tmpSpeciesvector
 */
acs_longInt environment::getTotalNumberOfMonomers()
{
    if(debugLevel == FINDERRORDURINGRUNTIME) cout << "\t\tenvironment::getTotalNumberOfMonomers start" << endl;

    acs_longInt tmpTotNumOfMols = 0;

    for(vector<species>::iterator tmpAllSpeciesIter = allSpecies.begin(); tmpAllSpeciesIter != allSpecies.end(); tmpAllSpeciesIter++)
    {
        if((tmpAllSpeciesIter->getAmount() > 0) && (tmpAllSpeciesIter->getSolubility() == SOLUBLE))
            tmpTotNumOfMols += tmpAllSpeciesIter->getAmount() * tmpAllSpeciesIter->getSequenceLength();
    }

    if(debugLevel == FINDERRORDURINGRUNTIME) cout << "\t\tenvironment::getTotalNumberOfMonomers end" << endl;

    return tmpTotNumOfMols;
}

/**
 Get the total number of charged molecules
*/
acs_int environment::getTotNumberOfChargedMols()
{
    if(debugLevel == FINDERRORDURINGRUNTIME) cout << "\t\tenvironment::getTotNumberOfChargedMold start" << endl;
    acs_int chargedMols = 0;
    for(vector<species>::iterator tmpAllSpeciesIter = allSpecies.begin(); tmpAllSpeciesIter != allSpecies.end(); tmpAllSpeciesIter++)
    {
            chargedMols += tmpAllSpeciesIter->getChargeMols();
    }
    if(debugLevel == FINDERRORDURINGRUNTIME) cout << "\t\tenvironment::getTotNumberOfChargedMold end" << endl;
    return chargedMols;
}

/**
 This Function populates two vectors containing all the nutrients and all the probabilities for the influx selected
 species from the firing disk according to the max length of the influx
 */
//TR void environment::getNutrientsFromTheFiringDisk()
//{
//	if(debugLevel == FINDERRORDURINGRUNTIME) cout << "\t\tenvironment::getNutrientsFromTheFiringDisk start" << endl;
//	
//	acs_double tmpProbability;
//
//	for(acs_int singleSpecies = 0; singleSpecies < (acs_int)firingDisk.size(); singleSpecies++)
//	{
//		// IF the species (contained in the firing disk) is longer (or at least equal) to the influx 
//		if(firingDisk.at(singleSpecies).getSequenceLength() <= influx)
//		{
//			nutrientsForInflux.push_back(firingDisk.at(singleSpecies).getID());
//			tmpProbability = createInitialConcentration(alphabet.length(), firingDisk.at(singleSpecies).getSequence().length()) / overallConcentration;
//			nutrientsProb2BeSelected.push_back(tmpProbability);
//		}
//	}
//	
//	if(debugLevel == FINDERRORDURINGRUNTIME) cout << "\t\tenvironment::getNutrientsFromTheFiringDisk end" << endl;
//}

/**
 This function print to monitor the content of the vectors nutrientsForInflux and nutrientsProb2BeSelected
 @version 1.0
 @date 2010-05-17
*/
void environment::printNutrientsAndProbability()
{
    if(debugLevel == FINDERRORDURINGRUNTIME) cout << "\t\tenvironment::printNutrientsAndProbability start" << endl;

    for(vector<acs_int>::iterator tmpAllNutrientsIter = nutrientsForInflux.begin(); tmpAllNutrientsIter != nutrientsForInflux.end(); tmpAllNutrientsIter++)
    {
        cout << *tmpAllNutrientsIter << "\t";
    }

    cout << endl;
	
    for(vector<acs_double>::iterator tmpAllNutrientsIter = nutrientsProb2BeSelected.begin(); tmpAllNutrientsIter != nutrientsProb2BeSelected.end(); tmpAllNutrientsIter++)
    {
        cout << *tmpAllNutrientsIter << "\t";
    }

    cout << endl;

    if(debugLevel == FINDERRORDURINGRUNTIME) cout << "\t\tenvironment::printNutrientsAndProbability end" << endl;
}

/**
 Nutrients amount fixing process. The amount of nutrients has to be fixed according to the initial theoretical distribution
 @version 1.0
 */
void environment::nutrientsAmountsFixing()
{
	if(debugLevel == FINDERRORDURINGRUNTIME) cout << "\t\t\tenvironment::nutrientsAmountsFixing start" << endl;
	
	// Cumulative temporary probability
	acs_double tmpTotProb = 0;
	
    for(vector<acs_double>::iterator tmpAllNutrientsIter = nutrientsProb2BeSelected.begin(); tmpAllNutrientsIter != nutrientsProb2BeSelected.end(); tmpAllNutrientsIter++)
	{
		tmpTotProb += *tmpAllNutrientsIter;
	}
	
    vector<acs_double>nutrientsProb2BeSelectedCOPY = nutrientsProb2BeSelected;
	//Clear wrong structure
	nutrientsProb2BeSelected.clear();
	
    for(vector<acs_double>::iterator tmpAllNutrientsIter = nutrientsProb2BeSelectedCOPY.begin(); tmpAllNutrientsIter != nutrientsProb2BeSelectedCOPY.end(); tmpAllNutrientsIter++)
	{
		nutrientsProb2BeSelected.push_back(*tmpAllNutrientsIter/tmpTotProb);
	}	

	if(debugLevel == FINDERRORDURINGRUNTIME) cout << "\t\t\tenvironment::nutrientsAmountsFixing end" << endl;
}


/**
 Perform Refill according to the time interval and the total amount of refill
 @version 2.4
 @date 2010.06.10
 @param tmpTimeSinceTheLastInFlux time elapsed since the last influx of at least one molcolule
 @param tmpMinimalTimeForOneMols time necessary to feed to the system one molecule
 @param tmp__RndDoubleGen random number generator
 */
bool environment::performRefill(acs_double tmpTimeSinceTheLastInFlux, acs_double tmpMinimalTimeForOneMols, MTRand& tmp__RndDoubleGen)
{
	if(debugLevel == FINDERRORDURINGRUNTIME) cout << "\tenvironment::performRefill start" << endl;
	bool refillFlag = true;
	
//TR	//IF THE MAX LENGTH TO REFILL IS GREATER THAN THE MAX LENGHT OF THE FIRING DISK
//	if(influx > initialMaxLength)
//	{
//		ExitWithError("performRefill", "ERROR: The max refill lenght is longer than the initial max length");
//		refillFlag = false;
//	}
	acs_int numberOfMolsToRefill = acsround(tmpTimeSinceTheLastInFlux/tmpMinimalTimeForOneMols);
	if(debugLevel >= SMALL_DEBUG ||  numberOfMolsToRefill > 1000)
	{
		cout << "\t\t\t\t|- Number of molecules to refill: " << numberOfMolsToRefill << " DT: " << tmpTimeSinceTheLastInFlux << " T: " << tmpMinimalTimeForOneMols << endl;
	}
	
	if(numberOfMolsToRefill > 0)
	{
            // SELECT NUTRIENTS FROM THE FIRING DISK
            for(acs_longInt singleRefill = 0; singleRefill < numberOfMolsToRefill; singleRefill++)
            {
                acs_longInt nutrientForInflux_ID = returnSelectionIdFromAWeightProbVectorAlreadyNormalized(nutrientsProb2BeSelected,tmp__RndDoubleGen);
                if(debugLevel >= HIGH_DEBUG)
                {
                        cout << "\t\t\t\t";
                        for(acs_int sngProb = 0; sngProb < (acs_int)nutrientsProb2BeSelected.size(); sngProb++)
                        {
                           cout << nutrientsProb2BeSelected.at(sngProb) << " ";
                        }
                        cout << endl;
                        cout << "\t\t\t\t|- Nutrient To refill: " << allSpecies.at(nutrientsForInflux.at(nutrientForInflux_ID)).getSequence() << "#" <<  nutrientForInflux_ID << endl;
                }

                if(nutrientForInflux_ID >= (unsigned)nutrientsForInflux.size())
                        cout << nutrientsForInflux.at(nutrientForInflux_ID) << " greater than the possible nutrients selection " << nutrientsForInflux.size()-1 << endl;

                acs_longInt IDspecies;
                try{
                    IDspecies = allSpecies.at(nutrientsForInflux.at(nutrientForInflux_ID)).getID();
                  }
                catch(exception&e){
                    cout << "allSpecies.at(nutrientsForInflux.at(nutrientForInflux_ID)).getID();" << endl;
                    cout << "Vectorsize " << allSpecies.size()<<" - position " << nutrientsForInflux.at(nutrientForInflux_ID) << endl;
                    cerr << "exceptioncaught:" << e.what() << endl;

                }

                if(debugLevel >= SMALL_DEBUG)
                {
                    cout << "\t\t\t\t|- " << singleRefill << "/" << numberOfMolsToRefill << " - Tot Number of Moles before refilling: " << getTotalNumberOfMolecules();
                    cout << " - Species to Refill: (" << allSpecies.at(IDspecies).getSequence() << "#" << IDspecies << ") - Tot Amount: " <<
                                    allSpecies.at(IDspecies).getAmount() << endl;
                }

                // ACTUAL INFLUX EVENT ----------------- increment specific species and the total number of species
                if(IDspecies >= (unsigned)allSpecies.size())
                    cout << IDspecies << " maggiore del massimo consentito che è " << allSpecies.size()-1 << endl;
                try{
                allSpecies.at(IDspecies).increment(volume);
                }
                catch(exception&e){
                    cout << " allSpecies.at(IDspecies).increment(volume);" << endl;
                    cout << "Vectorsize " << allSpecies.size()<<" - position " << IDspecies<< endl;
                    cerr << "exceptioncaught:" << e.what() << endl;

                }

                incMolProcedure(IDspecies);
                incSpeciesProcedure(IDspecies);

                if(debugLevel >= SMALL_DEBUG)
                {
                        cout << " - after refilling: " << getTotalNumberOfMolecules() << " - Tot Amount: " << allSpecies.at(IDspecies).getAmount() << endl;
                }
            }
	}else{
		refillFlag = false;
	}
	
	if(debugLevel == FINDERRORDURINGRUNTIME) cout << "\tenvironment::performRefill end" << endl;
	return refillFlag;
}

/**
 This function perform the deterministic dissociation of the complex
 @version 1.0
 @date 2013-09-18
 @param tmpTimeInterval time elapsed since the last reaction
 @param MTrand &tmp_RndDoubleGen random number generator
 */
bool environment::performDETComplexDissociation(acs_double tmpTimeInterval, MTRand& tmp_RndDoubleGen)
{
	if(debugLevel == FINDERRORDURINGRUNTIME) cout << "\tenvironment::performDETComplexDissociation start" << endl;
	bool cpxFlug = true;
	vector<acs_double> tempCpx;
	vector<acs_longInt> tempID;
	acs_double cpxDissRate = 0;
	acs_double totCpxDissRate = 0;
	acs_int cpxIntegerPart;
	acs_longInt selectedID;
	acs_int selectedPosition;
	acs_int traceID;
	for(vector<species>::iterator speciesIter = allSpecies.begin(); speciesIter < allSpecies.end(); speciesIter++)
	{
		// Insert complexes into a temporary list
		if((speciesIter->getComplexCutPnt() > 0) && (speciesIter->getAmount() > 0))
		{
			cpxDissRate = speciesIter->getAmount() * speciesIter->getComplexDegEnh() * tmpTimeInterval;
			totCpxDissRate += cpxDissRate;
			tempCpx.push_back(totCpxDissRate);
			tempID.push_back(speciesIter->getID());
			decimalComplexesToDissociate += cpxDissRate;
		}
	}
	if(decimalComplexesToDissociate > 1)
	{
		// Computer number of complexes to dissociate
			if(debugLevel == COMPLEXSTUFF) cout << "-----------" << endl;
			if(debugLevel == COMPLEXSTUFF) cout << decimalComplexesToDissociate << endl;


		// If the number of complexes to dissociate is greater than the actual number of complexes (due to round reason)
		// the cpxIntegerPart is set to the total number of complex molecules.
		if(decimalComplexesToDissociate >= numberOfCpxMols)
			decimalComplexesToDissociate = numberOfCpxMols;

		cpxIntegerPart = (acs_int)decimalComplexesToDissociate;

			if(debugLevel == COMPLEXSTUFF) cout << cpxIntegerPart << endl;

		decimalComplexesToDissociate = decimalComplexesToDissociate - cpxIntegerPart;

			if(debugLevel == COMPLEXSTUFF) cout << decimalComplexesToDissociate << endl;
			if(debugLevel == COMPLEXSTUFF) cout << tempCpx.size() << endl;

			if(decimalComplexesToDissociate < 0) decimalComplexesToDissociate = 0;
		for(acs_int i = 0; i < cpxIntegerPart; i++)
		{
			try{
				if(tempCpx.size() > 0)
				{
					// According to the propensity list the complex is randomly selected
					if(tempCpx.size() > 1)
					{
						selectedPosition = returnSelectionIdFromAWeightProbVector(tempCpx,totCpxDissRate,tmp_RndDoubleGen, __LINE__);
					}else{
						selectedPosition = 0;
					}
					selectedID = tempID.at(selectedPosition);
					if(debugLevel == COMPLEXSTUFF)
					{	cout << "-- REMOVING " << i+1 << " OF " << cpxIntegerPart
							<< "Pos " << selectedPosition << " - Sel " << selectedID << " - #" << allSpecies.at(selectedID).getAmount()
							<< " -cat "<< allSpecies.at(selectedID).getCatalyst_ID()
							<< " -sub " << allSpecies.at(selectedID).getSubstrate_ID() << endl;
					}
					// Perform complex dissociation
					if(!performComplexDissociation(selectedID,allSpecies.at(selectedID).getCatalyst_ID(),allSpecies.at(selectedID).getSubstrate_ID(),tmp_RndDoubleGen))
						ExitWithError("performComplexDissociation", "Problems during the COMPLEX DETERMINISTIC DISSOCIATION process");

					if(debugLevel == COMPLEXSTUFF) cout << "Residues - #" << allSpecies.at(selectedID).getAmount() << endl;

					// The list of possible dissociation is now recalculeted according to the previous dissociation
					traceID = 0; // -1 because list starts from 0
					vector<acs_double>::iterator backIterTempCpx = tempCpx.end()-1;
					while(backIterTempCpx >= tempCpx.begin()+selectedPosition)
					{
						*backIterTempCpx -= allSpecies.at(selectedID).getComplexDegEnh() * tmpTimeInterval;
						totCpxDissRate -= allSpecies.at(selectedID).getComplexDegEnh() * tmpTimeInterval;
						backIterTempCpx--;
						traceID++;
					}
					if(allSpecies.at(selectedID).getAmount() == 0)
					{
						if(debugLevel == COMPLEXSTUFF)
						{
							cout << "Processing the elementID " << tempCpx.size()-traceID  << ": " << *backIterTempCpx << endl;
							for(acs_int z = 0; z < tempCpx.size(); z++){cout << " - Pos " << z << " " << tempCpx.at(z) << " " << tempID.at(z) << endl;}
						}
						tempCpx.erase(tempCpx.end()-traceID);
						tempID.erase(tempID.end()-traceID);
						if(debugLevel == COMPLEXSTUFF)
						{
							cout << endl;
							for(acs_int z = 0; z < tempCpx.size(); z++){cout << " - Pos " << z << " " << tempCpx.at(z) << " " << tempID.at(z) << endl;}
						}
					}
				}else{
					cout << "\t\tNumber of complexes to dissociate: " << cpxIntegerPart
						 << " - Removing the " << i+1 << "th complex"
						 << " - Number of actual Complexes: " << tempCpx.size() << endl;
					decimalComplexesToDissociate = 0;
					//break;
					ExitWithError("performDETComplexDissociation","ERROR!!! Although the system compute possibile dissociation, noone complex is available");
				}

			}
				catch(exception&e)
			{
				cout << "ERROR in Deterministic Complex Dissociation, line " << __LINE__ << endl;
				cerr << "exceptioncaught:" << e.what() << endl;
				ExitWithError("performDETComplexDissociation","exceptionerrorthrown");
			}

		}
	}


	if(debugLevel == FINDERRORDURINGRUNTIME) cout << "\tenvironment::performDETComplexDissociation end" << endl;
	return cpxFlug;
}

/**
  This function perform the pseudo-deterministic molecules charging process
 @version 2.6
 @date 2011-02-24
 @param tmpTimeInterval time elapsed since the last reaction
 @param MTrand &tmp_RndDoubleGen random number generator
 */
bool environment::performDETMoleculesCharging(acs_double tmpTimeInterval, MTRand &tmp_RndDoubleGen)
{
    if(debugLevel == FINDERRORDURINGRUNTIME) cout << "\tenvironment::performDETMoleculesCharging start" << endl;
    bool effFlug = true;
    // ENERGY -------------------------------------
    if(energy == ENERGYBASED) // If energy is turned on molecules are charged
    {
       // cout << "******************************* BEGIN" << endl;
       // cout << "-- PRE" << endl;

        setNotChargedAndChargedSpeciesIDsAndAmounts(); // create two lists containing charged and unchaged molecules

		// If there are molecules to load
        if((!unchargedIDlist.empty()) && (!cumUnchargedAmountList.empty()))
        {

            acs_int nrgIntegerPart;
						
			// DIFFERENT FUNCTIONS TO LOAD DETERMINISTICALY HAVE BEEN IMPLMENTED 
            //acs_double actualAmountOfMolToLoad = (cumUnchargedAmountList.back() * K_nrg * tmpTimeInterval) + decimalMoleculesToLoad;
            //acs_double actualAmountOfMolToLoad =  ((ECConcentration * tmpTimeInterval * AVO) * (cumUnchargedAmountList.back()) * (K_nrg / (AVO * volume))) + decimalMoleculesToLoad;
			//acs_double actualAmountOfMolToLoad =  ((ECConcentration * AVO) * tmpTimeInterval) + decimalMoleculesToLoad;
			acs_double tempInEnergyMols = ECConcentration * AVO * tmpTimeInterval;
			acs_double tempInEnergyConc = tempInEnergyMols/(AVO*volume);
            acs_double notChargeMolsConc = cumUnchargedAmountList.back()/(AVO*volume);
			int loadedEnergyResidue = 0;
			
            acs_double actualAmountOfMolToLoad = (((K_nrg * tmpTimeInterval) * tempInEnergyConc * notChargeMolsConc) * (AVO * volume)) + decimalMoleculesToLoad;
			
            //TR cout << "................." << endl;
            //cout << "notChargeMolsConc" << notChargeMolsConc << endl;
            //cout << "tempInEnergyConc" << tempInEnergyConc << endl;
            //cout << "tmpTimeInterval" << tmpTimeInterval << endl;
            //cout << "knrg" << K_nrg << endl;
            //cout << "mol to load " << actualAmountOfMolToLoad << endl;
																	
            if((isinf(actualAmountOfMolToLoad)) || (actualAmountOfMolToLoad >= cumUnchargedAmountList.back()))
            { 
                nrgIntegerPart = cumUnchargedAmountList.back();
                loadedEnergyResidue = tempInEnergyMols - cumUnchargedAmountList.back();
				if(loadedEnergyResidue < 0)
				{
					decimalMoleculesToLoad += 0;
				}else{
					decimalMoleculesToLoad += loadedEnergyResidue;
				}
            }else{
				nrgIntegerPart = (acs_int)actualAmountOfMolToLoad;
				decimalMoleculesToLoad = actualAmountOfMolToLoad - (acs_double)nrgIntegerPart;
            }

            //cout << "\t2. after: "  << nrgIntegerPart << " " << decimalMoleculesToLoad << endl;

            // Declare temporary variable indicating the species to load
            acs_int tmpSpeciesToLoad;
            if(nrgIntegerPart > 0)
            {
                acs_int j;
                for(acs_int tmpSpeciesCnt = 0; tmpSpeciesCnt < nrgIntegerPart; tmpSpeciesCnt++)
                {
                    if((!unchargedIDlist.empty()) && (!cumUnchargedAmountList.empty()))
                    {
                        if(cumUnchargedAmountList.size() > 1)
                        {
                            tmpSpeciesToLoad = returnSelectionIdFromAWeightProbVector(cumUnchargedAmountList, cumUnchargedAmountList.back(), tmp_RndDoubleGen, __LINE__);
                        }else{
                            tmpSpeciesToLoad = 0;
                        }

                        // Molecule Charcing

                        if(!allSpecies.at(unchargedIDlist.at(tmpSpeciesToLoad)).chargeMol())
                        {
                           ExitWithError("PerformMoleculeEfflux", "Error during the charging process");
                        }else{
                           incOverallLoadedMolsCounter();
                        }

                        // When a species goes out the cumulative list is updated;
                        if(unchargedIDlist.size() > 0)
                        {
                            if(unchargedIDlist.size() == 1)
                            {
                                cumUnchargedAmountList.at(0)--;
                                if(cumUnchargedAmountList.at(0) == 0)
                                {
                                    cumUnchargedAmountList.clear();
                                    unchargedIDlist.clear();
                                }
                            }else{
                                j = tmpSpeciesToLoad;
                                while(j < (acs_int)unchargedIDlist.size())
                                {
                                    cumUnchargedAmountList.at(j)--;
                                    if(j == 0)
                                    {
                                        if(cumUnchargedAmountList.at(j) == 0)
                                        {
                                           //unchargedIDlist.removeAt(j);
                                           unchargedIDlist.erase(unchargedIDlist.begin()+j);
                                           //cumUnchargedAmountList.removeAt(j);
                                           cumUnchargedAmountList.erase(cumUnchargedAmountList.begin()+j);
                                        }else{
                                            j++;
                                        }
                                    }else{
                                       if(cumUnchargedAmountList.at(j) == cumUnchargedAmountList.at(j-1))
                                       {
                                          //unchargedIDlist.removeAt(j);
                                          unchargedIDlist.erase(unchargedIDlist.begin()+j);
                                          //cumUnchargedAmountList.removeAt(j);
                                          cumUnchargedAmountList.erase(cumUnchargedAmountList.begin()+j);
                                       }else{
                                          j++;
                                       }
                                    } // if(j == 0)
                                } // end while(j < (acs_int)unchargedIDlist.size())
                            } // if(unchargedIDlist.size() == 1)
                        } // end if(unchargedIDlist.size() > 0)
                        //cout << "post vector unchargedIDlist sizes " << unchargedIDlist.size() << " " << cumUnchargedAmountList.size() << endl;
                    } // end if((!unchargedIDlist.empty()) && (!cumUnchargedAmountList.empty()))
                } //  for(acs_int tmpSpeciesCnt = 0; tmpSpeciesCnt < nrgIntegerPart; tmpSpeciesCnt++)
            } // end if(nrgIntegerPart > 0)
        }else{
			decimalMoleculesToLoad += ECConcentration * tmpTimeInterval; // If all the molecules are charged all the avalaible energy carries wait for the next reaction
        } // end if((!unchargedIDlist.empty()) && (!cumUnchargedAmountList.empty()))
		
		// ---------------------------------|
		// Molecules Energy DECAY Procedure |
		// ---------------------------------|
        if((!chargedIDlist.empty()) && (!cumChargedAmountList.empty()))
		{
			acs_int nrgDecayIntegerPart;
                        acs_double actualAmountOfMolToUNLOAD = (K_nrg_decay * cumChargedAmountList.back() * tmpTimeInterval) + decimalMoleculesToUNLOAD;
			if(actualAmountOfMolToUNLOAD > 0)
            { 
				nrgDecayIntegerPart = (acs_int)actualAmountOfMolToUNLOAD;
				decimalMoleculesToUNLOAD = actualAmountOfMolToUNLOAD - (acs_double)nrgDecayIntegerPart;
            }else {
				nrgDecayIntegerPart = 0;
				decimalMoleculesToUNLOAD = 0;
			}

			// Declare temporary variable indicating the species to UNload
            acs_int tmpSpeciesToUNLoad;
			if(nrgDecayIntegerPart > 0)
			{
				acs_int z;
                for(acs_int tmpSpeciesCnt = 0; tmpSpeciesCnt < nrgDecayIntegerPart; tmpSpeciesCnt++)
                {
                    if((!chargedIDlist.empty()) && (!cumChargedAmountList.empty()))
                    {
                        if(cumChargedAmountList.size() > 1)
                        {
                            tmpSpeciesToUNLoad = returnSelectionIdFromAWeightProbVector(cumChargedAmountList, cumChargedAmountList.back(), tmp_RndDoubleGen, __LINE__);
                        }else{
                            tmpSpeciesToUNLoad = 0;
                        }
						
                        // Molecule Charcing
                        try{
                            if(!allSpecies.at(chargedIDlist.at(tmpSpeciesToUNLoad)).unchargeMol())
                            {
                                ExitWithError("PerformMoleculeEfflux", "Error during the decay process");
                            }else{
                                decOverallLoadedMolsCounter();
                            }
                        }
                        catch(exception&e)
                        {
                            cout << "if(!allSpecies.at(chargedIDlist.at(tmpSpeciesToUNLoad)).unchargeMol())" << endl;
                            cout << "Vectorsize " << allSpecies.size()<<" - position " << chargedIDlist.at(tmpSpeciesToUNLoad) << endl;
                            cerr << "exceptioncaught:" << e.what() << endl;
                            ExitWithError(" Molecule Charcing","exceptionerrorthrown");
                        }
                        // When a species goes out the cumulative list is updated;
                        if(chargedIDlist.size() > 0)
                        {
                            if(chargedIDlist.size() == 1)
                            {
                                cumChargedAmountList.at(0)--;
                                if(cumChargedAmountList.at(0) == 0)
                                {
                                    cumChargedAmountList.clear();
                                    chargedIDlist.clear();
                                }
                            }else{
                                z = tmpSpeciesToUNLoad;
                                while(z < (acs_int)chargedIDlist.size())
                                {
                                    cumChargedAmountList.at(z)--;
                                    if(z == 0)
                                    {
                                        if(cumChargedAmountList.at(z) == 0)
                                        {
                                            //chargedIDlist.removeAt(z);
                                            chargedIDlist.erase(chargedIDlist.begin()+z);
                                            //cumChargedAmountList.removeAt(z);
                                            cumChargedAmountList.erase(cumChargedAmountList.begin()+z);
                                        }else{
                                            z++;
                                        }
                                    }else{
                                            if(cumChargedAmountList.at(z) == cumChargedAmountList.at(z-1))
                                            {
                                                    //chargedIDlist.removeAt(z);
                                                    chargedIDlist.erase(chargedIDlist.begin()+z);
                                                    //cumChargedAmountList.removeAt(z);
                                                    cumChargedAmountList.erase(cumChargedAmountList.begin()+z);
                                            }else{
                                                    z++;
                                            }
                                    } // if(z == 0)
                                } // end while(z < (acs_int)chargedIDlist.size())
                            } // if(chargedIDlist.size() == 1)
                        } // end if(chargedIDlist.size() > 0)
                        //cout << "post vector chargedIDlist sizes " << chargedIDlist.size() << " " << cumChargedAmountList.size() << endl;
                    } // end if((!chargedIDlist.empty()) && (!cumChargedAmountList.empty()))
                } //  for(acs_int tmpSpeciesCnt = 0; tmpSpeciesCnt < nrgDecayIntegerPart; tmpSpeciesCnt++)
			}
		}else{ 
			decimalMoleculesToUNLOAD = 0;
        } // end if((!chargedIDlist.empty()) && (!cumChargedAmountList.empty()))
		
		
		
    } // end if(energy == ENERGYBASED)
    // END ENERGY ------------------------------

    if(debugLevel == FINDERRORDURINGRUNTIME) cout << "\tenvironment::performDETMoleculesCharging end" << endl;
    return effFlug;
}

/**
 This function populates to lists, the first containing all the living species ID and the
 second one containing the cumulative number of living species
 @version 2.5.2
 @date 2010-11-11
 */
void environment::setLivingSpeciesIDsAndAmounts()
{
    if(debugLevel == FINDERRORDURINGRUNTIME) cout << "\tenvironment::setLivingSpeciesIDsAndAmounts start" << endl;

    livingSpeciesIDsList.clear(); // Clean the livingSpeciesIDs list
    cumLivingAmountsList.clear(); // Clean living Species Amount
    acs_longInt cumulativeNumberOfSpecies = 0; // Initialise the cumulative number of species
    bool toBeInserted = true;
	
	// cout << ".................." << endl;

    if(!allSpecies.empty())
    {
        for(vector<species>::iterator tmpAllSpeciesIter = allSpecies.begin(); tmpAllSpeciesIter != allSpecies.end(); tmpAllSpeciesIter++)
        {
            toBeInserted = true;
            if(maxLOut > 0){if(tmpAllSpeciesIter->getSequenceLength() <= maxLOut) toBeInserted = true; else toBeInserted = false;}
			
			if(toBeInserted)
			{
				if((tmpAllSpeciesIter->getAmount() > 0) && (tmpAllSpeciesIter->getSolubility() == SOLUBLE))
				{
					// cout << tmpAllSpeciesIter->getID() << " " << tmpAllSpeciesIter->getSequence() << " " << tmpAllSpeciesIter->getAmount() <<  endl;
					livingSpeciesIDsList.push_back(tmpAllSpeciesIter->getID());
					cumulativeNumberOfSpecies += tmpAllSpeciesIter->getAmount();
					cumLivingAmountsList.push_back((acs_double)cumulativeNumberOfSpecies);              
				}
			}
        }
    }


    if(debugLevel == FINDERRORDURINGRUNTIME) cout << "\tenvironment::setLivingSpeciesIDsAndAmounts end" << endl;
}

/**
 This function populates to lists, the first one contains the IDs of the uncharged mols and the second one
 contains the cumulative amount of uncharged molecules
 @version 2.5.3
 @date 2011-02-22
 */
void environment::setNotChargedAndChargedSpeciesIDsAndAmounts()
{
    if(debugLevel == FINDERRORDURINGRUNTIME) cout << "\tenvironment::setNotChargedAndChargedSpeciesIDsAndAmounts start" << endl;

    unchargedIDlist.clear(); // Clean Living Species Uncharged ID list
    cumUnchargedAmountList.clear(); // Clean Living specie cumulative uncharged List
    acs_longInt cumulativeNumberOfUnchargedSpecies = 0; // Initialise the cumulative number of species
	
	// Clear and declare the same variables as above, but for the charged species. 
	chargedIDlist.clear();
	cumChargedAmountList.clear();
	acs_longInt cumulativeNumberOfChargedSpecies = 0; 

    if(!allSpecies.empty())
    {
        for(vector<species>::iterator tmpAllSpeciesIter = allSpecies.begin(); tmpAllSpeciesIter != allSpecies.end(); tmpAllSpeciesIter++)
        {
            if((tmpAllSpeciesIter->getAmount() > 0) && (tmpAllSpeciesIter->getSolubility() == SOLUBLE) && (tmpAllSpeciesIter->getComplexCutPnt() == 0))
            {
                if(energy == ENERGYBASED)
                {
                    if((tmpAllSpeciesIter->getNOTchargeMols() > 0) && (tmpAllSpeciesIter->getEnergizable()))
                    {
                        unchargedIDlist.push_back(tmpAllSpeciesIter->getID());
                        cumulativeNumberOfUnchargedSpecies += tmpAllSpeciesIter->getNOTchargeMols();
                        cumUnchargedAmountList.push_back((acs_double)cumulativeNumberOfUnchargedSpecies);
                    }
					if((tmpAllSpeciesIter->getChargeMols() > 0) && (tmpAllSpeciesIter->getEnergizable()))
					{
						chargedIDlist.push_back(tmpAllSpeciesIter->getID());
						cumulativeNumberOfChargedSpecies += tmpAllSpeciesIter->getChargeMols();
						cumChargedAmountList.push_back((acs_double)cumulativeNumberOfChargedSpecies);
					}
                }
            }
        }
    }

    if(debugLevel == FINDERRORDURINGRUNTIME) cout << "\tenvironment::setNotChargedAndChargedSpeciesIDsAndAmounts end" << endl;
}

/**
 Perform vector unchargedIDlist update removing a new charge molecule
 vectors unchargedIDlist and cumUnchargedAmountList are involved
 @version 1.0
 @date 2010-10-10
 @param acs_int tmpSpeciesID Specie to uncharge
 */
bool environment::removeChargeMolFromList(acs_int tmpSpeciesID)
{
    if(debugLevel == FINDERRORDURINGRUNTIME) cout << "\tenvironment::removeChargeMolFromList start" << endl;

    bool chFlag = false;
    acs_int counter = 0;
    bool fromHere = false;
    for(vector<acs_int>::iterator tmpAllSpeciesIter = unchargedIDlist.begin(); tmpAllSpeciesIter != unchargedIDlist.end(); tmpAllSpeciesIter++)
    {
        if(*tmpAllSpeciesIter == tmpSpeciesID)
            fromHere = true;

        if(fromHere)
        {
            try{
            cumUnchargedAmountList.at(counter)++;
            }
            catch(exception&e)
            {
                cout << "cumUnchargedAmountList.at(counter)++;" << endl;
                cout << "Vectorsize " << cumUnchargedAmountList.size()<<" - position " << counter << endl;
                cerr << "exceptioncaught:" << e.what() << endl;
                ExitWithError("cumUnchargedAmountList","exceptionerrorthrown");
            }
            chFlag = true;
        }

        counter++;
    }
    // If the species is not present it is added at the end of the list with the last value of the list + 1
    if(!fromHere)
    {
        unchargedIDlist.push_back(tmpSpeciesID);
        cumUnchargedAmountList.push_back(cumUnchargedAmountList.back()+1);
        chFlag = true;
    }

    if(debugLevel == FINDERRORDURINGRUNTIME) cout << "\tenvironment::removeChargeMolFromList end" << endl;
    return chFlag;
}

/**
 Perform vector unchargedIDlist update adding a new charge molecule
 vector unchargedIDlist and cumUnchargedAmountList are involved
 @version 1.0
 @date 2010-10-10
 @param acs_int tmpSpeciesID Specie to charge
 */
bool environment::addChargeMolToList(acs_int tmpSpeciesID)
{
    if(debugLevel == FINDERRORDURINGRUNTIME) cout << "\tenvironment:addChargeMolToList start" << endl;

    bool chFlag = false;
    acs_int counter = 0;
    bool fromHere = false;
    bool flagRemove = false;
    for(vector<acs_int>::iterator tmpAllSpeciesIter = unchargedIDlist.begin(); tmpAllSpeciesIter != unchargedIDlist.end(); tmpAllSpeciesIter++)
    {
        try{
        if(*tmpAllSpeciesIter == tmpSpeciesID)
        {
            fromHere = true;
            if(cumUnchargedAmountList.at(counter) == 1)
                flagRemove = true;
        }

        if(fromHere)
        {
            cumUnchargedAmountList.at(counter)--;
            if(!flagRemove)
                chFlag = true;
        }
        }
        catch(exception&e)
        {
        cout<<" cumUnchargedAmountList.at(counter)--" << endl;
        cout << "Vectorsize " <<cumUnchargedAmountList.size()<<" - position " <<counter<< endl;
        cerr << "exceptioncaught:" << e.what() << endl;
        ExitWithError("cumUnchargedAmountList","exceptionerrorthrown");
        }

        counter++;
    }
    // If the species is not present it is added at the end of the list with the last value of the list + 1
    if(flagRemove)
    {
        //unchargedIDlist.removeAt(tmpSpeciesID);
        unchargedIDlist.erase(unchargedIDlist.begin()+tmpSpeciesID);
        //cumUnchargedAmountList.removeAt(tmpSpeciesID);
        cumUnchargedAmountList.erase(cumUnchargedAmountList.begin()+tmpSpeciesID);
        chFlag = true;
    }

    if(debugLevel == FINDERRORDURINGRUNTIME) cout << "\tenvironment::addChargeMolToList end" << endl;
    return chFlag;
}

/**
 This function perform the deterministic molecules efflux process
 @version 2.5
 @date 2010-06-25
 @param tmpTimeInterval time elapsed since the last reaction
 @param MTrand &tmp_RndDoubleGen random number generator
 */
bool environment::performMoleculesEfflux(acs_double tmpTimeInterval, MTRand &tmp_RndDoubleGen)
{
    if(debugLevel == FINDERRORDURINGRUNTIME) cout << "\tenvironment::performMoleculesEfflux start" << endl;
    bool effFlug = true;
    if(debugLevel == SMALL_DEBUG){cout << "\t\t\t\t|- Number of molecules: BEFORE efflux: " << getTotalNumberOfMolecules();}
	
    // Create the list of all the living species IDs and Amounts and ENERGY as well
    setLivingSpeciesIDsAndAmounts();
	
	// EFFLUX
	
	// energy efflux
	// If the system is open some energy molecules go away => moleculeDecay_KineticConstant
	if(nrgBoolFlag == ENERGYBASED)
	{
		acs_double tempNrgOut = decimalMoleculesToLoad * moleculeDecay_KineticConstant * tmpTimeInterval;
		if(tempNrgOut > 0)
		{
			if((decimalMoleculesToLoad - tempNrgOut) > 0)
			{
				decimalMoleculesToLoad = decimalMoleculesToLoad - tempNrgOut;
			}
			else{
				decimalMoleculesToLoad = 0;
			}
		}
		
		// ... and some EC molecules decay --> K_nrg_decay
		tempNrgOut = decimalMoleculesToLoad * K_nrg_decay * tmpTimeInterval;
		if(tempNrgOut > 0)
		{
			if((decimalMoleculesToLoad - tempNrgOut) > 0)
			{
				decimalMoleculesToLoad = decimalMoleculesToLoad - tempNrgOut;
			}
			else{
				decimalMoleculesToLoad = 0;
			}		
		}
	}
	
	// If there are molecules...
    if((!cumLivingAmountsList.empty()) && (!livingSpeciesIDsList.empty()))
    {
		
        // Compute the number of species to efflux according to the delta T elapsed from the last reaction
        acs_int integerPart;		
		
        acs_double actualAmountOfMoltoEfflux = (cumLivingAmountsList.back() * moleculeDecay_KineticConstant * tmpTimeInterval) + decimalMoleculesToEfflux;
		
        // cout << cumLivingAmountsList.back() << " * " << moleculeDecay_KineticConstant<< " * " << tmpTimeInterval << " + " << decimalMoleculesToEfflux << " = " << actualAmountOfMoltoEfflux << endl;
		
        if((isinf(actualAmountOfMoltoEfflux)) || (actualAmountOfMoltoEfflux > cumLivingAmountsList.back()))
        {
            integerPart = cumLivingAmountsList.back();
			decimalMoleculesToEfflux += 0;
        }else{
			integerPart = (acs_int)actualAmountOfMoltoEfflux;
			decimalMoleculesToEfflux = actualAmountOfMoltoEfflux - (acs_double)integerPart;
        }
		
        acs_longInt speciesToEfflux;
		
        // For each molecule to efflux
        if(integerPart > 0)
        {
            acs_int i;
            for(acs_longInt singleEfflux = 0; singleEfflux < integerPart; singleEfflux++)
            {
                if((!cumLivingAmountsList.empty()) && (!livingSpeciesIDsList.empty()))
                {
                    // Select species ID (ID of the livingSpeciesIDsList list) to efflux
                    if(livingSpeciesIDsList.size() > 1)
                    {
                        speciesToEfflux = returnSelectionIdFromAWeightProbVector(cumLivingAmountsList, cumLivingAmountsList.back(),tmp_RndDoubleGen, __LINE__);
                    }else{
                        speciesToEfflux = 0;
                    }
                    // Molecule efflux will be livingSpeciesIDsList at the position speciesToEfflux
                    if(debugLevel == SMALL_DEBUG)
                    {
						cout << "\t\t\tDT\t" << tmpTimeInterval << "\tSpecies to efflux\t" << livingSpeciesIDsList.at(speciesToEfflux) <<
						" \tMols <- efflux: " << allSpecies.at(livingSpeciesIDsList.at(speciesToEfflux)).getAmount() << endl;
                    }
                    try{
                    // ACTUAL EFFLUX PROCEDURE
                    if(!performMoleculeEfflux(livingSpeciesIDsList.at(speciesToEfflux), tmp_RndDoubleGen))
						ExitWithError("performMoleculeEfflux", "Problems during the EFFLUX process");
                    }
                    catch(exception&e)
                    {
                    cout<<"!performMoleculeEfflux(livingSpeciesIDsList.at(speciesToEfflux), tmp_RndDoubleGen))" << endl;
                    cout << "Vectorsize " <<livingSpeciesIDsList.size()<<" - position " <<speciesToEfflux<< endl;
                    cerr << "exceptioncaught:" << e.what() << endl;
                    ExitWithError("speciesToEfflux","exceptionerrorthrown");
                    }

                    if(livingSpeciesIDsList.size() > 0)
                    {
                        if(livingSpeciesIDsList.size() == 1)
                        {
                            cumLivingAmountsList.at(0)--;
                            if(cumLivingAmountsList.at(0) == 0)
                            {
                                cumLivingAmountsList.clear();
                                livingSpeciesIDsList.clear();
                            }
                        }else{
                            i = speciesToEfflux;
                            while(i < (acs_int)livingSpeciesIDsList.size())
                            {
                                cumLivingAmountsList.at(i)--;
                                if(i == 0)
                                {
                                    if(cumLivingAmountsList.at(i) == 0)
                                    {
                                        //livingSpeciesIDsList.removeAt(i);
                                        livingSpeciesIDsList.erase(livingSpeciesIDsList.begin()+i);
                                        //cumLivingAmountsList.removeAt(i);
                                        cumLivingAmountsList.erase(cumLivingAmountsList.begin()+i);
                                    }else{
                                        i++;
                                    }
                                }else{
                                    if(cumLivingAmountsList.at(i) == cumLivingAmountsList.at(i-1))
                                    {
                                        //livingSpeciesIDsList.removeAt(i);
                                        livingSpeciesIDsList.erase(livingSpeciesIDsList.begin()+i);
                                        //cumLivingAmountsList.removeAt(i);
                                        cumLivingAmountsList.erase(cumLivingAmountsList.begin()+i);
                                    }else{
                                        i++;
                                    }
                                } // if(i == 0)
                            } // end while(i < (acs_int)livingSpeciesIDsList.size())
                        }
                    }
					
                    if(debugLevel == SMALL_DEBUG)
                    {
						cout << "\tMols -> efflux: " << allSpecies.at(livingSpeciesIDsList.at(speciesToEfflux)).getAmount() << endl;
                    }
                }else{
                    ExitWithError("performMoleculesEfflux", "No molecules to efflux");
                } //  if((!cumLivingAmountsList.empty()) && (!livingSpeciesIDsList.empty()))
            } //  for(acs_longInt singleEfflux = 0; singleEfflux < integerPart; singleEfflux++)
        } // END if(integerPart > 0)
    } // if((!cumLivingAmountsList.empty()) && (!livingSpeciesIDsList.empty()))
	
    if(debugLevel == SMALL_DEBUG){cout << " - AFTER efflux: " << getTotalNumberOfMolecules() << endl;}
    if(debugLevel == FINDERRORDURINGRUNTIME) cout << "\tenvironment::performMoleculesEfflux end" << endl;
    return effFlug;
}


/**
 Perform MOLECULE EFFLUX reaction
 @version 2.5.1
 @date 2010-06-27
 @param acs_longInt tmpSpecies Species ID
 @param MTRand& tmp__RndDoubleGen random generator
 */
bool environment::performMoleculeEfflux(acs_longInt tmpSpecies, MTRand& tmp__RndDoubleGen)
{
	if(debugLevel == FINDERRORDURINGRUNTIME) cout << "\t\tenvironment::performMoleculeEfflux start" << endl;
	bool reactionFlag = false;
	
	if(debugLevel >= SMALL_DEBUG)
	{
		cout << "\t\t\t|- Molecule to Efflux: " << allSpecies.at(tmpSpecies).getSequence() << " - Before efflux process [#" << allSpecies.at(tmpSpecies).getAmount() << "]";
	}
	
	// Complex (loaded or not loaded) to dissociate is chosen according to the total amount of loaded and not loaded complexes present in the system
    acs_double ratioLoadedSpecies;

    try{
    ratioLoadedSpecies= allSpecies.at(tmpSpecies).getChargeMols() / allSpecies.at(tmpSpecies).getAmount();
    }
    catch(exception&e)
    {
    cout<<" ratioLoadedSpecies= allSpecies.at(tmpSpecies).getChargeMols() / allSpecies.at(tmpSpecies).getAmount()" << endl;
    cout << "Vectorsize " <<allSpecies.size()<<" - tmpSpecies " << tmpSpecies << endl;
    cerr << "exceptioncaught:" << e.what() << endl;
    ExitWithError("ratioLoadedSpecies","exceptionerrorthrown");
    }

    acs_double rndNumb = tmp__RndDoubleGen();

	if(rndNumb > ratioLoadedSpecies) // rndNumb has to be stricktly greater 
	{
		// DECREMENT TOTAL AMOUNT OF THIS species or complex (not loaded)
                if(allSpecies.at(tmpSpecies).getNOTchargeMols() < 1)
		{
			ExitWithError("performMoleculeEfflux", "There are not avalaible not loaded species to efflux");
		}else{
                        allSpecies.at(tmpSpecies).decrement(volume);
			if(allSpecies.at(tmpSpecies).getComplexCutPnt() == 0) // is the species is not a complex
			{
				// Here the concentrationFixed species proprieties is not evalute since it works only with "close" systems
				decMolSpeciesProcedure(tmpSpecies); // Change overall amount of species and molecules
			}else{
				decCpxProcedure(tmpSpecies);
			}
			reactionFlag = true;
		}
	}else{
		// DECREMENT TOTAL AMOUNT OF THIS species or complex (loaded)
                if(allSpecies.at(tmpSpecies).getChargeMols() < 1)
		{
			ExitWithError("performMoleculeEfflux", "There are not avalaible loaded species to efflux");
		}else{
			//cout << "------------------------" << endl;
			//cout << allSpecies.at(tmpSpecies).getChargeMols() << " of " << allSpecies.at(tmpSpecies).getAmount() << " || ";
			unchargeMolProcess(tmpSpecies);
                        allSpecies.at(tmpSpecies).decrement(volume);
			//cout << " " << allSpecies.at(tmpSpecies).getChargeMols() << " of " << allSpecies.at(tmpSpecies).getAmount() << endl;
			if(allSpecies.at(tmpSpecies).getComplexCutPnt() == 0) // is the species is not a complex
			{
				decMolSpeciesProcedure(tmpSpecies); // Change overall amount of species and molecules
			}else{
				decCpxProcedure(tmpSpecies);
			}
			reactionFlag = true;
		}
	}
	
	if(debugLevel >= SMALL_DEBUG)
	{
		cout << " - after efflux process [#" << allSpecies.at(tmpSpecies).getAmount() << "]" << endl;
	}
	
	if(debugLevel == FINDERRORDURINGRUNTIME) cout << "\t\tenvironment::performMoleculeEfflux end" << endl;
	return reactionFlag;	
}

/**
 Perform ENERGY EFFLUX reaction
 @version 2.4.1
 @date 2010-06-27
 @param MTRand& tmp__RndDoubleGen random number generator
 */

bool environment::performEnergyEfflux(MTRand& tmp__RndDoubleGen)
{
	if(debugLevel == FINDERRORDURINGRUNTIME) cout << "\t\tenvironment::performEnergyEfflux start" << endl;
	bool reactionFlag = false;
	
	if(debugLevel >= SMALL_DEBUG)
	{
		cout << "\t\t\t|- CARRIER to Efflux: " << nrgCarrier.getSequence() << " - Before efflux process [#" << nrgCarrier.getAmount() << "]";
	}
	
	// Complex (loaded or not loaded) to dissociate is chosen according to the total amount of loaded and not loaded complexes present in the system
	acs_double ratioLoadedCarriers = nrgCarrier.getChargeMols() / nrgCarrier.getAmount();
	acs_double rndNumb = tmp__RndDoubleGen();
	if(rndNumb > ratioLoadedCarriers) // If an unloaded energy carrier has to be throw away
	{
		
		// DECREMENT TOTAL AMOUNT OF THIS species (not loaded)
		if(nrgCarrier.getNOTchargeMols() < 1)
		{
			ExitWithError("performEnergyEfflux", "There are not avalaible not loaded carriers to efflux");
		}else{
			nrgCarrier.decrement(volume);
			reactionFlag = true;
		}
	}else{
		// DECREMENT TOTAL AMOUNT OF THIS species (loaded)
		if(nrgCarrier.getChargeMols() < 1)
		{
			ExitWithError("performEnergyEfflux", "There are not avalaible loaded carriers to efflux");
		}else{
			nrgCarrier.unchargeMol();
			nrgCarrier.decrement(volume);
			reactionFlag = true;
		}
	}
	
	if(debugLevel >= SMALL_DEBUG)
	{
		cout << " - after carrier efflux process [#" << nrgCarrier.getAmount() << "]" << endl;
	}
	
	if(debugLevel == FINDERRORDURINGRUNTIME) cout << "\t\tenvironment::performEnergyEfflux end" << endl;
	return reactionFlag;
	
}

/**
 Perform the reaction after the Gillespie computation
 @version 1.1
 @date 2013-10-29
 @param acs_longInt reaction_u reaction ID in Gillespie structure
 @param MTRand& tmp_RndDoubleGen Random number generator
 @param acs_int tmp_ActGEN Current generation
 @param acs_int tmp_ActSIM Current generation 
 @param acs_int tmp_ActSTEP Current step (reaction)
 @param string Storing directory
 */
bool environment::performReaction(acs_longInt reaction_u, MTRand& tmp_RndDoubleGen, acs_int tmp_ActGEN, 
								  acs_int tmp_ActSIM, acs_int tmp_ActSTEP, string tmp_StoringPath)
{
	if(debugLevel == FINDERRORDURINGRUNTIME) cout << "\tenvironment::performReaction start" << endl;

	bool reactionFlag = true;
	
	// ^-^-^-^-^-^-^-^-^-^-^-^-^-^-^-^-^-^-^-^-^-^-^-^-^-^-^-^-^-^-^-^-^-^-^-
	// PERFORM OPERATIONS ACCORDING TO THE REACTION
	// ^-^-^-^-^-^-^-^-^-^-^-^-^-^-^-^-^-^-^-^-^-^-^-^-^-^-^-^-^-^-^-^-^-^-^-
    try{
        allGillespieScores.at(reaction_u);
    }
    catch(exception&e)
    {
    cout<<" allGillespieScores.at(reaction_u)" << endl;
    cout << "Vectorsize " <<allGillespieScores.size()<<" - position " << reaction_u << endl;
    cerr << "exceptioncaught:" << e.what() << endl;
    if(debugLevel >= SMALL_DEBUG) printGillespieStructure();
    ExitWithError("environment::performReaction","exceptionerrorthrown");
    }

	switch (allGillespieScores.at(reaction_u).getIdReactionType()) 
    {
		case CONDENSATION:
			if(!performCondensation(allGillespieScores.at(reaction_u).getMolI(), //catalyst
									allGillespieScores.at(reaction_u).getMolII(), // substrate 
									allGillespieScores.at(reaction_u).getMolIII(), // product
									allGillespieScores.at(reaction_u).getMolIV(), // complex	
									allGillespieScores.at(reaction_u).getReactionID(),
									allGillespieScores.at(reaction_u).getCatalysisID(), 
									tmp_RndDoubleGen))
			{
			if(debugLevel >= RUNNING_VERSION)
				cout << "\t\t\t|- !*!*!* There's no molecule availability to perform CONDENSATION reaction..." << endl;
			}	
			// STORE ON FILE TIMES EVERY fileTimesSaveInterval seconds (if at least something happen)
			if((getActualTime() > (getFileTimesSavingInterval() + internalTimesStoredCounter)) ||
					(getActualTime() == 0) || (getFileTimesSavingInterval() == 0))
			{
				saveReactionsParametersSTD(tmp_ActGEN, tmp_ActSIM, tmp_ActSTEP, tmp_StoringPath, allGillespieScores.at(reaction_u).getIdReactionType(),
										allGillespieScores.at(reaction_u).getMolI(), allGillespieScores.at(reaction_u).getMolIII(),
										allSpecies.at(allGillespieScores.at(reaction_u).getMolIV()).getSubstrate_ID(), allGillespieScores.at(reaction_u).getMolII());

				saveLivingSpeciesIDSTD(tmp_ActGEN, tmp_ActSIM, tmp_ActSTEP, tmp_StoringPath);
				saveLivingSpeciesAmountSTD(tmp_ActGEN, tmp_ActSIM, tmp_StoringPath);
				saveLivingSpeciesConcentrationSTD(tmp_ActGEN, tmp_ActSIM, tmp_StoringPath);
				if(getActualTime() > 0)
				{
					internalTimesStoredCounter = internalTimesStoredCounter + getFileTimesSavingInterval();
				}
			}
			break;
		case ENDO_CONDENSATION:
			if(!perform_endo_Condensation(allGillespieScores.at(reaction_u).getMolI(), //catalyst
									allGillespieScores.at(reaction_u).getMolII(), // substrate
									allGillespieScores.at(reaction_u).getMolIII(), // product
									allGillespieScores.at(reaction_u).getMolIV(), // complex
									allGillespieScores.at(reaction_u).getNRGside(), // Energy species side, complex, substrate or both
									allGillespieScores.at(reaction_u).getReactionID(),
									allGillespieScores.at(reaction_u).getCatalysisID(), 
									tmp_RndDoubleGen))
			{
				if(debugLevel >= RUNNING_VERSION)
					cout << "\t\t\t|- !*!*!* There's no molecule availability to perform ENDO CONDENSATION reaction..." << endl;
            }
			if((getActualTime() > (getFileTimesSavingInterval() + internalTimesStoredCounter)) ||
			   (getActualTime() == 0))
			{
				try{
				saveReactionsParametersSTD(tmp_ActGEN, tmp_ActSIM, tmp_ActSTEP, tmp_StoringPath, allGillespieScores.at(reaction_u).getIdReactionType(),
										allGillespieScores.at(reaction_u).getMolI(), allGillespieScores.at(reaction_u).getMolIII(),
										allSpecies.at(allGillespieScores.at(reaction_u).getMolIV()).getSubstrate_ID(), allGillespieScores.at(reaction_u).getMolII());
				}
				catch(exception&e)
				{
				cout<<" saveReactionsParametersSTD(tmp_ActGEN, tmp_ActSIM, tmp_ActSTEP, tmp_StoringPath, allGillespieScores.at(reaction_u).getIdReactionType(allGillespieScores.at(reaction_u).getMolI(), allGillespieScores.at(reaction_u).getMolIII(),allSpecies.at(allGillespieScores.at(reaction_u).getMolIV()).getSubstrate_ID(), allGillespieScores.at(reaction_u).getMolII());" << endl;
				cout << __LINE__ << endl;
				cerr << "exceptioncaught:" << e.what() << endl;
				ExitWithError("","exceptionerrorthrown");
				}

				saveLivingSpeciesIDSTD(tmp_ActGEN, tmp_ActSIM, tmp_ActSTEP, tmp_StoringPath);
				saveLivingSpeciesAmountSTD(tmp_ActGEN, tmp_ActSIM, tmp_StoringPath);
				saveLivingSpeciesConcentrationSTD(tmp_ActGEN, tmp_ActSIM, tmp_StoringPath);
				if(getActualTime() > 0)
				{
					internalTimesStoredCounter = internalTimesStoredCounter + getFileTimesSavingInterval();
				}
			}
			break;
		case CLEAVAGE:
			if(!performCleavage(allGillespieScores.at(reaction_u).getMolI(), // Substrate
								allGillespieScores.at(reaction_u).getMolII(), // Product I
								allGillespieScores.at(reaction_u).getMolIII(), // product II
								allGillespieScores.at(reaction_u).getReactionID(), // reaction ID 
								allGillespieScores.at(reaction_u).getCatalysisID(), // Catalysis ID 
								tmp_RndDoubleGen))
				ExitWithError("performCleavage", "Problems during the CLEAVAGE reaction");	
			
			if((getActualTime() > (getFileTimesSavingInterval() + internalTimesStoredCounter)) || (getActualTime() == 0))
			{
				saveReactionsParametersSTD(tmp_ActGEN, tmp_ActSIM, tmp_ActSTEP, tmp_StoringPath, allGillespieScores.at(reaction_u).getIdReactionType(),
										   allGillespieScores.at(reaction_u).getMolIV(), allGillespieScores.at(reaction_u).getMolI(), allGillespieScores.at(reaction_u).getMolII(),
										   allGillespieScores.at(reaction_u).getMolIII());

				saveLivingSpeciesIDSTD(tmp_ActGEN, tmp_ActSIM, tmp_ActSTEP, tmp_StoringPath);
				saveLivingSpeciesAmountSTD(tmp_ActGEN, tmp_ActSIM, tmp_StoringPath);
				saveLivingSpeciesConcentrationSTD(tmp_ActGEN, tmp_ActSIM, tmp_StoringPath);
				if(getActualTime() > 0){ internalTimesStoredCounter = internalTimesStoredCounter + getFileTimesSavingInterval(); }
			}
			break;
		case ENDO_CLEAVAGE:
			if(!perform_endo_Cleavage(allGillespieScores.at(reaction_u).getMolI(), // Substrate
									  allGillespieScores.at(reaction_u).getMolII(),	// Product One
									  allGillespieScores.at(reaction_u).getMolIII(), // Product Two
									  allGillespieScores.at(reaction_u).getNRGside(), // Energy species side, complex, substrate or both
									  allGillespieScores.at(reaction_u).getReactionID(), // reaction ID
									  allGillespieScores.at(reaction_u).getCatalysisID(), // Catalysis ID
									  tmp_RndDoubleGen))
                                ExitWithError("perform_endo_Cleavage", "Problems during the ENDOCLEAVAGE reaction");
			
			if((getActualTime() > (getFileTimesSavingInterval() + internalTimesStoredCounter)) || (getActualTime() == 0))
			{
				saveReactionsParametersSTD(tmp_ActGEN, tmp_ActSIM, tmp_ActSTEP, tmp_StoringPath, allGillespieScores.at(reaction_u).getIdReactionType(),
										   allGillespieScores.at(reaction_u).getMolIV(), allGillespieScores.at(reaction_u).getMolI(), allGillespieScores.at(reaction_u).getMolII(),
										   allGillespieScores.at(reaction_u).getMolIII());

				saveLivingSpeciesIDSTD(tmp_ActGEN, tmp_ActSIM, tmp_ActSTEP, tmp_StoringPath);
				saveLivingSpeciesAmountSTD(tmp_ActGEN, tmp_ActSIM, tmp_StoringPath);
				saveLivingSpeciesConcentrationSTD(tmp_ActGEN, tmp_ActSIM, tmp_StoringPath);
				if(getActualTime() > 0) { internalTimesStoredCounter = internalTimesStoredCounter + getFileTimesSavingInterval(); }
			}
			break;
		case COMPLEXFORMATION:
			if(!performComplexFormation(allGillespieScores.at(reaction_u).getMolI(), // Catalyst
                                        allGillespieScores.at(reaction_u).getMolII(), // Substrate
                                        allGillespieScores.at(reaction_u).getMolIII(), // catalysis ID
                                        allGillespieScores.at(reaction_u).getMolIV(), // Second Substrate
										tmp_RndDoubleGen))
			{
				if(debugLevel >= RUNNING_VERSION)
					cout << "\t\t\t|- !*!*!* There's no molecule availability to perform COMPLEX FORMATION reaction..." << endl;
				
				ExitWithError("performComplexFormation", "Problems during the COMPLEX COMPLETION process");
			}			
			break;
		case ENDO_COMPLEXFORMATION:
			if(!perform_endo_ComplexFormation(allGillespieScores.at(reaction_u).getMolI(), // Catalyst
											  allGillespieScores.at(reaction_u).getMolII(), // Substrate
											  allGillespieScores.at(reaction_u).getMolIII(), // catalysis ID
											  allGillespieScores.at(reaction_u).getMolIV(), // Second Substrate
                                              allGillespieScores.at(reaction_u).getNRGside(), // Energy species side, complex, substrate or both
											  tmp_RndDoubleGen))
			{
				if(debugLevel >= RUNNING_VERSION)
					cout << "\t\t\t|- !*!*!* There's no molecule availability to perform ENDO COMPLEX FORMATION reaction..." << endl;
				ExitWithError("performComplexFormation", "Problems during the ENDO COMPLEX COMPLETION process");
			}			
			break;
		case COMPLEXDEGRADATION:
			if(!performComplexDissociation(allGillespieScores.at(reaction_u).getMolI(), 
										   allGillespieScores.at(reaction_u).getMolII(), 
										   allGillespieScores.at(reaction_u).getMolIII(),
										   tmp_RndDoubleGen))
				ExitWithError("performComplexDissociation", "Problems during the COMPLEX DISSOCIATION process");				
			break;
			
		case SPECIESDECAY:
			if(!performMoleculeEfflux(allGillespieScores.at(reaction_u).getMolI(), tmp_RndDoubleGen))
				ExitWithError("performMoleculeEfflux", "Problems during the EFFLUX process");			
			break;
			
		case SPONTANEOUS_CLEAVAGE:
			if(!performSpontaneousCleavage( allGillespieScores.at(reaction_u).getReactionID(), tmp_RndDoubleGen))
				ExitWithError("performSpontaneousCleavage", "Problem with spontaneous cleavage reaction");
			if((getActualTime() > (getFileTimesSavingInterval() + internalTimesStoredCounter)) || (getActualTime() == 0))
			{
				saveReactionsParametersSTD(tmp_ActGEN, tmp_ActSIM, tmp_ActSTEP, tmp_StoringPath, allGillespieScores.at(reaction_u).getIdReactionType(),
										   allGillespieScores.at(reaction_u).getMolIV(), allGillespieScores.at(reaction_u).getMolI(), allGillespieScores.at(reaction_u).getMolII(),
										   allGillespieScores.at(reaction_u).getMolIII());

				saveLivingSpeciesIDSTD(tmp_ActGEN, tmp_ActSIM, tmp_ActSTEP, tmp_StoringPath);
				saveLivingSpeciesAmountSTD(tmp_ActGEN, tmp_ActSIM, tmp_StoringPath);
				saveLivingSpeciesConcentrationSTD(tmp_ActGEN, tmp_ActSIM, tmp_StoringPath);
				if(getActualTime() > 0) {internalTimesStoredCounter = internalTimesStoredCounter + getFileTimesSavingInterval();}
			}
			break;
			
		case SPONTANEOUS_CONDENSATION:
			if(!performSpontaneousCondensation( allGillespieScores.at(reaction_u).getReactionID(), tmp_RndDoubleGen))
				ExitWithError("performSpontaneousCondensation", "Problem with spontaneous condensation reaction");
			if((getActualTime() > (getFileTimesSavingInterval() + internalTimesStoredCounter)) || (getActualTime() == 0) || (getFileTimesSavingInterval() == 0))
			{
				saveReactionsParametersSTD(tmp_ActGEN, tmp_ActSIM, tmp_ActSTEP, tmp_StoringPath, allGillespieScores.at(reaction_u).getIdReactionType(),
										allGillespieScores.at(reaction_u).getMolI(), allGillespieScores.at(reaction_u).getMolIII(),
										allSpecies.at(allGillespieScores.at(reaction_u).getMolIV()).getSubstrate_ID(), allGillespieScores.at(reaction_u).getMolII());

				saveLivingSpeciesIDSTD(tmp_ActGEN, tmp_ActSIM, tmp_ActSTEP, tmp_StoringPath);
				saveLivingSpeciesAmountSTD(tmp_ActGEN, tmp_ActSIM, tmp_StoringPath);
				saveLivingSpeciesConcentrationSTD(tmp_ActGEN, tmp_ActSIM, tmp_StoringPath);
				if(getActualTime() > 0){ internalTimesStoredCounter = internalTimesStoredCounter + getFileTimesSavingInterval(); }
			}
			break;
			
		case ENERGYEFFLUX:
			if (!performEnergyEfflux(tmp_RndDoubleGen))
				ExitWithError("performEnergyEfflux", "Problems during the ENERGY EFFLUX process");
			break;

        //TR case PHOSPHORILATION:
        //	if (!performPhosphorilation(allGillespieScores.at(reaction_u).getMolI()))
        //		ExitWithError("performPhosphorilation", "Problems during the PHOSPHORILATION process");
        //	break;
	}
	
	if(debugLevel == FINDERRORDURINGRUNTIME) cout << "\tenvironment::performReaction end" << endl;
	
	return reactionFlag;
}

/**
 Perform CONDENSATION reaction
 @version 1.2
 @date 2011-02-12
 @param acs_longInt tmpComplex Complex ID
 @param acs_longInt tmpSubstrate Substrate ID
 @param acs_longInt tmpProduct Product ID
 @param MTRand& tmp__RndDoubleGen random generator
 */
bool environment::performCondensation(acs_longInt tmpCatalyst, acs_longInt tmpSubstrate, acs_longInt tmpProduct, 
                                      acs_longInt tmpComplex, acs_longInt tmpIdReaction, acs_longInt tmpIdCatalysis,
                                      MTRand& tmp__RndDoubleGen)
{
	if(debugLevel == FINDERRORDURINGRUNTIME) cout << "\t\tenvironment::performCondensation start" << endl;

	bool reactionFlag = false;
	acs_longInt tmpSubAmount;
	acs_longInt tmpCpxAmount;

    try{
	if(nrgBoolFlag == ENERGYBASED)
	{
		tmpSubAmount = allSpecies.at(tmpSubstrate).getNOTchargeMols();
		tmpCpxAmount = allSpecies.at(tmpComplex).getNOTchargeMols();
	}else {
		tmpSubAmount = allSpecies.at(tmpSubstrate).getAmount();
		tmpCpxAmount = allSpecies.at(tmpComplex).getAmount();		
	}

	}
    catch(exception&e)
    {
    cout<<"tmpSubAmount = allSpecies.at(tmpSubstrate).getNOTchargeMols(); tmpCpxAmount = allSpecies.at(tmpComplex).getNOTchargeMols(); " << endl;
    cerr << "exceptioncaught:" << e.what() << endl;
    ExitWithError("","exceptionerrorthrown");
    }


	// check whether all the species amounts to decrement are greater than 1
	if(tmpSubAmount == 0)
		ExitWithError("performCondensation", "Substrate not Avalaible!!!");
	if(tmpCpxAmount == 0)
		ExitWithError("performCondensation", "Complex not Avalaible!!!");		

	
		
	// the catalyst come to be free
	if (debugLevel >= SMALL_DEBUG) 
	{
		cout << "CONDENSATION BEFORE complexes: " << tmpCpxAmount
			 << " - substrates: " << tmpSubAmount
                         << " - catalyst: " << allSpecies.at(tmpCatalyst).getAmount() << endl;
	}
	
    // REACTION!!!
    try{
        // concentration fixed is checked within the function
        allSpecies.at(tmpCatalyst).increment(volume);
    }
    catch(exception&e)
    {
        cout<<" allSpecies.at(tmpCatalyst).increment(volume);" << endl;
        cout << "Vectorsize " <<allSpecies.size()<<" - position " << tmpCatalyst << endl;
        cerr << "exceptioncaught:" << e.what() << endl;
        ExitWithError("","exceptionerrorthrown");
    }

    if(!allSpecies.at(tmpCatalyst).getConcentrationFixed())
    {
        incMolProcedure(tmpCatalyst); // Increment total number of molecules and species
        incSpeciesProcedure(tmpCatalyst);
    }

	// The substrate and the complex are consumed (the first substrate has been already decremented during the complex formation reaction
	if ((tmpSubAmount > 0) && (tmpCpxAmount > 0))
	{
		// REACTION!!!
		
        allSpecies.at(tmpSubstrate).decrement(volume);
        if(!allSpecies.at(tmpSubstrate).getConcentrationFixed())
            decMolSpeciesProcedure(tmpSubstrate); // decrement total number of molecules and, if so, species
	
        allSpecies.at(tmpComplex).decrement(volume);
        if(!allSpecies.at(tmpComplex).getConcentrationFixed())
            decCpxProcedure(tmpComplex); // decrement total number of complexes token and, if so, types
		
		reactionFlag = true;
		if (debugLevel >= SMALL_DEBUG) 
		{
                        cout << "CONDENSATION AFTER complexes: " << allSpecies.at(tmpComplex).getNOTchargeMols()
                        << " - substrates: " << allSpecies.at(tmpSubstrate).getNOTchargeMols()
                        << " - catalyst: " << allSpecies.at(tmpCatalyst).getNOTchargeMols() << endl;
		}	
	}else{
		ExitWithError("performCondensation", "Complex or second substrate are not avalaible");
	}
	// Update catalysis and Reactions structure with a secure check on ID reaction
	if(allCatalysis.at(tmpIdCatalysis).getReactionID() == allReactions.at(tmpIdReaction).getID())
	{
                allCatalysis.at(tmpIdCatalysis).updateTotAmount();
                allReactions.at(tmpIdReaction).updateTotEvents();
	}else{
		ExitWithError("performCondensation", "Problem with reaction ID coherence");
	}
	// EVALUATE PRODUCT
	if(!newSpeciesEvaluationIII(tmpProduct, tmp__RndDoubleGen))
		ExitWithError("newSpeciesEvaluationIII", "Problems during the new species evalutation");
	if(debugLevel == FINDERRORDURINGRUNTIME) cout << "\t\tenvironment::performCondensation end" << endl;
	
	// Update reaction counter
	incCondensationCounter();
	return reactionFlag;
}

/**
 Perform ENDO CONDENSATION reaction
 @version 1.2
 @date 2011-02-12
 @param tmpCatalyst Catalyst (bound in the complex) ID
 @param tmpSubstrate second (or first, depend on the condensation type) substrate ID
 @param tmpProduct product ID
 @param acs_longInt tmpComplex Complex ID
 @param acs_longInt tmpSubstrate Substrate ID
 @param acs_longInt tmpProduct Product ID
 @param MTRand& tmp__RndDoubleGen random generator
 */
bool environment::perform_endo_Condensation(acs_longInt tmpCatalyst, acs_longInt tmpSubstrate, acs_longInt tmpProduct, 
											acs_longInt tmpComplex, acs_int tmpNRGside, acs_longInt tmpIdReaction, 
											acs_longInt tmpIdCatalysis, MTRand& tmp__RndDoubleGen)
{
	if(debugLevel == FINDERRORDURINGRUNTIME) cout << "\t\tenvironment::perform_endo_Condensation start" << endl;
	
	bool reactionFlag = false;
	acs_longInt tmpSubAmount;
	acs_longInt tmpCpxAmount;


	
	// According to the NRG side position the right amount of molecules belonging to the complexes and substrates is selected
	switch (tmpNRGside) {
		case COMPLEXLOAD:
			tmpSubAmount = allSpecies.at(tmpSubstrate).getNOTchargeMols();
			tmpCpxAmount = allSpecies.at(tmpComplex).getChargeMols();
			break;
		case SUBSTRATELOAD:
			tmpSubAmount = allSpecies.at(tmpSubstrate).getChargeMols();
			tmpCpxAmount = allSpecies.at(tmpComplex).getNOTchargeMols();
			break;			
		default:
			tmpSubAmount = allSpecies.at(tmpSubstrate).getChargeMols(); 
			tmpCpxAmount = 	allSpecies.at(tmpComplex).getChargeMols();		
			break;
	}
	
	// check whether all the species amount of not loaded molecules to decrement are greater than 1
	if(tmpSubAmount == 0)
		ExitWithError("performCondensation", "Not loaded Substrate not Avalaible!!!");
	if(tmpCpxAmount == 0)
		ExitWithError("performCondensation", "Loaded Complex not Avalaible!!!");
	
	if (debugLevel >= SMALL_DEBUG) 
	{
		cout << "ENDOCONDENSATION BEFORE complexes: " << tmpCpxAmount
		<< " - substrates: " << tmpSubAmount
                << " - catalyst: " << allSpecies.at(tmpCatalyst).getAmount() << endl;
	}
	
	// REACTION! the catalyst come to be free
        allSpecies.at(tmpCatalyst).increment(volume);
    if(!allSpecies.at(tmpCatalyst).getConcentrationFixed())
    {
        incMolProcedure(tmpCatalyst); // Increment total number of molecules and species
        incSpeciesProcedure(tmpCatalyst);
    }
			
	// The substrate and the complex are consumed (the first substrate has been already decremented during the complex formation reaction
	if((tmpSubAmount > 0) && (tmpCpxAmount > 0))
	{		
            // REACTION!!!
		// According to the NRG side position complexes and/or substrates are uncharged 
		switch (tmpNRGside) {
			case COMPLEXLOAD:
				unchargeMolProcess(tmpComplex);
				break;
			case SUBSTRATELOAD:
				unchargeMolProcess(tmpSubstrate);
				break;			
			default:
				unchargeMolProcess(tmpComplex); 
				unchargeMolProcess(tmpSubstrate);		
				break;
		}		
		
                try{
                    allSpecies.at(tmpSubstrate).decrement(volume);
                    if(!allSpecies.at(tmpSubstrate).getConcentrationFixed())
                        decMolSpeciesProcedure(tmpSubstrate); // decrement total number of molecules and, if so, species
                }catch(exception& e){
                    cout << "Error in perform_endo_condensation allSpecies.at(tmpSubstrate).decrement(volume)..." << endl;
                    cerr << "exception caught: " << e.what() << endl;
                    ExitWithError("perform_endo_condensation","exception error thrown");
                }

                try{
                    allSpecies.at(tmpComplex).decrement(volume);
                    if(!allSpecies.at(tmpComplex).getConcentrationFixed())
                        decCpxProcedure(tmpComplex); // decrement total number of complexes token and, if so, types
                }catch(exception& e){
                    cout << "Error in perform_endo_condensation allSpecies.at(tmpComplex).decrement(volume)..." << endl;
                    cerr << "exception caught: " << e.what() << endl;
                    ExitWithError("perform_endo_condensation","exception error thrown");
                }

		reactionFlag = true;
		if (debugLevel >= SMALL_DEBUG)
		{
                        cout << "ENDOCONDENSATION AFTER complexes: " << allSpecies.at(tmpComplex).getAmount()
                        << " - substrates: " << allSpecies.at(tmpSubstrate).getAmount()
                        << " - catalyst: " << allSpecies.at(tmpCatalyst).getAmount() << endl;
		}
	}else{
            ExitWithError("perform_endo_Condensation", "Complex or second substrate are not avalaible");
	}
	// Update catalysis and Reactions structure with a secure check on ID reaction
    if(allCatalysis.at(tmpIdCatalysis).getReactionID() == allReactions.at(tmpIdReaction).getID())
    {       try{
            allCatalysis.at(tmpIdCatalysis).updateTotAmount();
            allReactions.at(tmpIdReaction).updateTotEvents();
        }
        catch(exception&e)
        {
            cout<<" allCatalysis.at(tmpIdCatalysis).updateTotAmount(); and  allReactions.at(tmpIdReaction).updateTotEvents();" << endl;
            cout << "Vectorsize allCatalysis " <<allCatalysis.size()<<" - position " << tmpIdCatalysis << endl;
            cout << "Vectorsize allReactions " <<allReactions.size()<<" - position " << tmpIdReaction << endl;
            cerr << "exceptioncaught:" << e.what() << endl;
            ExitWithError("Update catalysis and Reactions structure with a secure check on ID reaction","exceptionerrorthrown");
        }
    }else{
            ExitWithError("performCondensation", "Problem with reaction ID coherence");
	}
	// EVALUATE PRODUCT
	if(!newSpeciesEvaluationIII(tmpProduct, tmp__RndDoubleGen))
            ExitWithError("newSpeciesEvaluationIII", "Problems during the new species evalutation");
	
	if(debugLevel == FINDERRORDURINGRUNTIME) cout << "\t\tenvironment::perform_endo_Condensation end" << endl;
	
	// Update reaction counter
	incEndoCondensationCounter();
	return reactionFlag;
}

/**
 Perform CLEAVAGE reaction
 @version 1.2
 @date 2011.02.12
 @param acs_longInt tmpSubstrate Substrate ID
 @param acs_longInt tmpProduct_I Product 1 ID
 @param acs_longInt tmpProduct_II Product 2 ID
 @param tmpIdReaction Rections ID
 @param tmpIdCatalysis Catalysis ID
 @param MTRand& tmp__RndDoubleGen random number generator
 */
bool environment::performCleavage(acs_longInt tmpSubstrate, acs_longInt tmpProduct_I, acs_longInt tmpProduct_II,
									acs_longInt tmpIdReaction, acs_longInt tmpIdCatalysis, MTRand& tmp__RndDoubleGen)
{
	if(debugLevel == FINDERRORDURINGRUNTIME) cout << "\t\tenvironment::performCleavage start" << endl;

	bool reactionFlag = false;
	bool avalability = true;
	
	if(debugLevel >= SMALL_DEBUG)
		cout << "\t\t\t|- Species " <<  tmpSubstrate << " will be cleaved" << endl;

    try{
	
	// substrate decrement
	if((!allSpecies.at(tmpSubstrate).getAmount() > 0) || (!allSpecies.at(allCatalysis.at(tmpIdCatalysis).getCat()).getAmount() > 0))
	{
		ExitWithError("performCleavage", "Substrate or catalyst not Avalaible!!!");
	}else{
		if(allCatalysis.at(tmpIdCatalysis).getCat() == tmpSubstrate)
		{
			if(allSpecies.at(tmpSubstrate).getAmount() < 2)
				avalability = false;
		}
		if(avalability)
		{
            allSpecies.at(tmpSubstrate).decrement(volume);
            if(!allSpecies.at(tmpSubstrate).getConcentrationFixed())
                decMolSpeciesProcedure(tmpSubstrate); // Update overall number of species and molecules
			reactionFlag = true;
		}else{
			ExitWithError("performCleavage", "Substrate not Avalaible, when catalyst and substrate are the same, at least of that should be avalable");
		}
	}

    }
    catch(exception&e)
    {
    cout<<"" << endl;
    cerr << "exceptioncaught:" << e.what() << endl;
    ExitWithError("something wrong in substrate decrement within performCleavag method","exceptionerrorthrown");
    }

    try{

        // Update catalysis and Reactions structure with a secure check on ID reaction
        if(allCatalysis.at(tmpIdCatalysis).getReactionID() == allReactions.at(tmpIdReaction).getID())
        {
            allCatalysis.at(tmpIdCatalysis).updateTotAmount();
            allReactions.at(tmpIdReaction).updateTotEvents();
        }else{
        	cout << "Catalysis ID: " << allCatalysis.at(tmpIdCatalysis).getCatId()
        		 << " - Catalysis Reaction ID: " << allCatalysis.at(tmpIdCatalysis).getReactionID()
        	     << " - Reaction ID: " << allReactions.at(tmpIdReaction).getID() << endl;
            ExitWithError("performCleavage", "Problem with reaction ID coherence");
        }

    }
    catch(exception&e)
    {
        cout<<"Update catalysis and Reactions structure with a secure check on ID reaction" << endl;
        cerr << "exceptioncaught:" << e.what() << endl;
        ExitWithError("something wrong in within performCleavag method","exceptionerrorthrown");
    }

    try{
        // EVALUATE PRODUCT 1
        if(!newSpeciesEvaluationIII(tmpProduct_I, tmp__RndDoubleGen))
            ExitWithError("newSpeciesEvaluationIII", "Problems during the new species evalutation");
    }
    catch(exception&e)
    {
        cout<<" if(!newSpeciesEvaluation(allSpecies.at(tmpProduct_I).getSequence(), tmp__RndDoubleGen))" << endl;
        cerr << "exceptioncaught:" << e.what() << endl;
        ExitWithError("something wrong in evaluating product 1 in performCleavage method","exceptionerrorthrown");
    }

    try{
        // EVALUATE PRODUCT 2
        if(!newSpeciesEvaluationIII(tmpProduct_II, tmp__RndDoubleGen))
            ExitWithError("newSpeciesEvaluationIII", "Problems during the new species evalutation");
    }
    catch(exception&e)
    {
        cout<<" if(!newSpeciesEvaluation(allSpecies.at(tmpProduct_II).getSequence(), tmp__RndDoubleGen))" << endl;
        cerr << "exceptioncaught:" << e.what() << endl;
        ExitWithError("something wrong in evaluating product 2 in performCleavage method","exceptionerrorthrown");
    }

	if(debugLevel == FINDERRORDURINGRUNTIME) cout << "\t\tenvironment::performCleavage end" << endl;

        // Update reaction counter
        incCleavageCounter();
	return reactionFlag;
}

/**
 Perform ENDO_CLEAVAGE reaction
 @version 1.2
 @date 2010.11.08
 @param acs_longInt tmpSubstrate Substrate ID
 @param acs_longInt tmpProduct_I Product 1 ID
 @param acs_longInt tmpProduct_II Product 2 ID
 @param acs_int tmpNrgTarget Energy Target --> catalyst loaded, substrate loaded or both
 @param tmpIdReaction Rections ID
 @param tmpIdCatalysis Catalysis ID
 @param MTRand& tmp__RndDoubleGen random number generator
 */
bool environment::perform_endo_Cleavage(acs_longInt tmpSubstrate, acs_longInt tmpProduct_I, acs_longInt tmpProduct_II,
										acs_int tmpNRGside, acs_longInt tmpIdReaction, acs_longInt tmpIdCatalysis, MTRand& tmp__RndDoubleGen)
{
    if(debugLevel == FINDERRORDURINGRUNTIME) cout << "\t\tenvironment::perform_endo_Cleavage start" << endl;

    bool reactionFlag = false;

    if(debugLevel >= SMALL_DEBUG)
            cout << "\t\t\t|- Species " <<  tmpSubstrate << " will be cleaved" << endl;

    try{

        switch (tmpNRGside) {
        case CATALYSTLOAD:
            if((!allSpecies.at(allCatalysis.at(tmpIdCatalysis).getCat()).getChargeMols() > 0) || (!allSpecies.at(tmpSubstrate).getNOTchargeMols() > 0))
            {
                ExitWithError("perform_endo_Cleavage", "Loaded substrate not Avalaible!!!");
            }else{
                unchargeMolProcess(allCatalysis.at(tmpIdCatalysis).getCat());
                allSpecies.at(tmpSubstrate).decrement(volume);
                if(!allSpecies.at(tmpSubstrate).getConcentrationFixed())
                    decMolSpeciesProcedure(tmpSubstrate); // Update overall number of species and molecules
                reactionFlag = true;
            }
            break;
        case SUBSTRATELOAD:
            if((!allSpecies.at(allCatalysis.at(tmpIdCatalysis).getCat()).getNOTchargeMols() > 0) || (!allSpecies.at(tmpSubstrate).getChargeMols() > 0))
            {
                ExitWithError("perform_endo_Cleavage", "Loaded substrate not Avalaible!!!");
            }else{
                unchargeMolProcess(tmpSubstrate);
                allSpecies.at(tmpSubstrate).decrement(volume);
                if(!allSpecies.at(tmpSubstrate).getConcentrationFixed())
                    decMolSpeciesProcedure(tmpSubstrate); // Update overall number of species and molecules
                reactionFlag = true;
            }
            break;
        default:
            if((!allSpecies.at(allCatalysis.at(tmpIdCatalysis).getCat()).getChargeMols() > 0) || (!allSpecies.at(tmpSubstrate).getChargeMols() > 0))
            {
                ExitWithError("perform_endo_Cleavage", "Loaded substrate not Avalaible!!!");
            }else{
                unchargeMolProcess(allCatalysis.at(tmpIdCatalysis).getCat());
                unchargeMolProcess(tmpSubstrate);
                allSpecies.at(tmpSubstrate).decrement(volume);
                if(!allSpecies.at(tmpSubstrate).getConcentrationFixed())
                    decMolSpeciesProcedure(tmpSubstrate); // Update overall number of species and molecules
                reactionFlag = true;
            }
            break;
        }

    }
    catch(exception&e)
    {
        cout<<"switch (tmpNRGside)" << endl;
        cerr << "exceptioncaught:" << e.what() << endl;
        ExitWithError("something wrong in switch block of method perform_endo_Cleavage","exceptionerrorthrown");
    }

    try{

        // Update catalysis and Reactions structure with a secure check on ID reaction
        if(allCatalysis.at(tmpIdCatalysis).getReactionID() == allReactions.at(tmpIdReaction).getID())
        {
            allCatalysis.at(tmpIdCatalysis).updateTotAmount();
            allReactions.at(tmpIdReaction).updateTotEvents();
        }else{
            ExitWithError("perform_endo_Cleavage", "Problem with reaction ID coherence");
        }

    }
    catch(exception&e)
    {
        cout<<"Update catalysis and Reactions structure with a secure check on ID reaction" << endl;
        cerr << "exceptioncaught:" << e.what() << endl;
        ExitWithError("something wrong in if block in perform_endo_Cleavage method","exceptionerrorthrown");
    }

    try{
        // EVALUATE PRODUCT 1
        if(!newSpeciesEvaluationIII(tmpProduct_I, tmp__RndDoubleGen))
            ExitWithError("newSpeciesEvaluationIII", "Problems during the new species evalutation");
    }
    catch(exception&e)
    {
        cout<<"if(!newSpeciesEvaluationIII(allSpecies.at(tmpProduct_I).getSequence(), tmp__RndDoubleGen))" << endl;
        cerr << "exceptioncaught:" << e.what() << endl;
        ExitWithError("something wrong in evaluating product 1 in  form_endo_Cleavage method","exceptionerrorthrown");
    }


    try{
        // EVALUATE PRODUCT 2
        if(!newSpeciesEvaluationIII(tmpProduct_II, tmp__RndDoubleGen))
            ExitWithError("newSpeciesEvaluationIII", "Problems during the new species evalutation");
    }
    catch(exception&e)
    {
        cout<<"if(!newSpeciesEvaluationIII(allSpecies.at(tmpProduct_II).getSequence(), tmp__RndDoubleGen))" << endl;
        cerr << "exceptioncaught:" << e.what() << endl;
        ExitWithError("something wrong in evaluating product 2 in  form_endo_Cleavage method","exceptionerrorthrown");
    }


    if(debugLevel == FINDERRORDURINGRUNTIME) cout << "\t\tenvironment::perform_endo_Cleavage end" << endl;

    // Update reaction counter
    incEndoCleavageCounter();
    return reactionFlag;
}

/**
 Perform COMPLEX FORMATION reaction
 @version 1.2
 @date 2011.02.13
 @param acs_longInt tmpCatalyst Catalyst ID
 @param acs_longInt tmpSubstrate Substrate ID
 @tmpCatalysisID Catalysis ID
 @param MTRand& tmp__RndDoubleGen random generator
 */
bool environment::performComplexFormation(acs_longInt tmpCatalyst, acs_longInt tmpSubstrate, acs_longInt tmpCatID,
										  acs_longInt tmpSecSub, MTRand& tmp__RndDoubleGen)
{
	if(debugLevel == FINDERRORDURINGRUNTIME) cout << "\t\tenvironment::performComplexFormation start" << endl;
	
	bool reactionFlag = false;
	bool avalability = true;
        string strErr = "";
	
	if(debugLevel >= SMALL_DEBUG)
            cout << "\t\t\t|- Complex between " <<  tmpCatalyst << "(" << allSpecies.at(tmpCatalyst).getSequence()  << "#" << allSpecies.at(tmpCatalyst).getNOTchargeMols() << ") and "
             << tmpSubstrate << "(" << allSpecies.at(tmpSubstrate).getSequence()  << "#" << allSpecies.at(tmpSubstrate).getNOTchargeMols() << ") will be formed -> CP: "
             << allSpecies.at(tmpCatalyst).getSequence().length() << endl;

    try{
        // IF CATALYST AND SUBSTRATE ARE EQUAL AT LEAST 2 MOLECULES HAVE TO BE PRESENT
        if(tmpCatalyst == tmpSubstrate)
        {
            if(allSpecies.at(tmpSubstrate).getNOTchargeMols() < 2)
            {
                avalability = false;
            }
        }
    }
    catch(exception&e)
    {
        cout<<"IF CATALYST AND SUBSTRATE ARE EQUAL AT LEAST 2 MOLECULES HAVE TO BE PRESENT" << endl;
        cerr << "exceptioncaught:" << e.what() << endl;
        ExitWithError("something wrong in if(tmpCatalyst == tmpSubstrate) block in performComplexFormation ","exceptionerrorthrown");
    }

	if(avalability)
	{
        string tmpComplex;

        try{
            // CREATE COMPLEX
            tmpComplex = allSpecies.at(tmpCatalyst).getSequence() + allSpecies.at(tmpSubstrate).getSequence();
        }
        catch(exception&e)
        {
            cout<<"tmpComplex = allSpecies.at(tmpCatalyst).getSequence() + allSpecies.at(tmpSubstrate).getSequence();" << endl;
            cerr << "exceptioncaught:" << e.what() << endl;
            ExitWithError("something wrong in creating complex in performComplexFormation","exceptionerrorthrown");
        }

        if(allSpecies.at(tmpCatalyst).getNOTchargeMols() < 1 || allSpecies.at(tmpSubstrate).getNOTchargeMols() < 1)
		{
				avalability = false;
		}

		bool tempNrgType = ENDOERGONIC;
		
		if(avalability)
		{
            allSpecies.at(tmpCatalyst).decrement(volume);
            if(!allSpecies.at(tmpCatalyst).getConcentrationFixed())
                decMolSpeciesProcedure(tmpCatalyst); // Update overall number of species and molecules

                        allSpecies.at(tmpSubstrate).decrement(volume);
            if(!allSpecies.at(tmpSubstrate).getConcentrationFixed())
                decMolSpeciesProcedure(tmpSubstrate); // Update overall number of species and molecules
			reactionFlag = true;
		}else{
			ExitWithError("performComplexFormation", strErr);
		}
		// EVALUATE COMPLEX (allSpecies.at(tmpCatalyst).getSequence().length() is definitly the cutting point)
		if(debugLevel == COMPLEXSTUFF){cout << "PERFORMCOMPLEXFORMATION |- Catalysis: " << tmpCatID << " sec sub: " << tmpSecSub << endl;}
		if(!complexEvaluation(tmpComplex, tmp__RndDoubleGen, allSpecies.at(tmpCatalyst).getSequence().length(),
                              tmpCatalyst, tmpSubstrate, tmpCatID, tmpSecSub, tempNrgType))
				ExitWithError("complexEvaluation", "Problems during the complex evalutation");
	}
	
	incCpxFormCounter();

	if(debugLevel == FINDERRORDURINGRUNTIME) cout << "\t\tenvironment::performComplexFormation end" << endl;
		
	return reactionFlag;	
}
/**
 Perform ENDO COMPLEX FORMATION reaction
 @version 1.2
 @date 2011.04.13
 @param acs_longInt tmpCatalyst Catalyst ID
 @param acs_longInt tmpSubstrate Substrate ID
 @tmpCatalysisID Catalysis ID
 @param MTRand& tmp__RndDoubleGen random generator
 */
bool environment::perform_endo_ComplexFormation(acs_longInt tmpCatalyst, acs_longInt tmpSubstrate, acs_longInt tmpCatID,
												acs_longInt tmpSecSub, acs_int tmpNRGSide, MTRand& tmp__RndDoubleGen)
{
	if(debugLevel == FINDERRORDURINGRUNTIME) cout << "\t\tenvironment::perform_endo_ComplexFormation start" << endl;
	
	bool reactionFlag = false; //This flag has to becames true to validate the reaction
        //string strErr = "";
	string tmpComplex;
	
	if(debugLevel >= SMALL_DEBUG)
		cout << "\t\t\t|- Complex between " <<  tmpCatalyst << "(" << allSpecies.at(tmpCatalyst).getSequence()  << "#" << allSpecies.at(tmpCatalyst).getNOTchargeMols() << ") and "
		<< tmpSubstrate << "(" << allSpecies.at(tmpSubstrate).getSequence()  << "#" << allSpecies.at(tmpSubstrate).getChargeMols() << ") will be formed -> CP: "
		<< allSpecies.at(tmpCatalyst).getSequence().length() << endl;

    try{
        // CREATE COMPLEX
        tmpComplex = allSpecies.at(tmpCatalyst).getSequence() + allSpecies.at(tmpSubstrate).getSequence();
    }
    catch(exception&e)
    {
        cout<<"tmpComplex = allSpecies.at(tmpCatalyst).getSequence() + allSpecies.at(tmpSubstrate).getSequence();" << endl;
        cerr << "exceptioncaught:" << e.what() << endl;
        ExitWithError("something wrong in creating complex in perform_endo_ComplexFormation","exceptionerrorthrown");
    }

        if(tmpNRGSide == SUBSTRATELOAD)
        {
            if((allSpecies.at(tmpCatalyst).getNOTchargeMols() > 0) && (allSpecies.at(tmpSubstrate).getChargeMols() > 0))
            {
				allSpecies.at(tmpCatalyst).decrement(volume);
				if(!allSpecies.at(tmpCatalyst).getConcentrationFixed())
					decMolSpeciesProcedure(tmpCatalyst); // Update overall number of species and molecules

				unchargeMolProcess(tmpSubstrate);
				allSpecies.at(tmpSubstrate).decrement(volume);
				if(!allSpecies.at(tmpSubstrate).getConcentrationFixed())
					decMolSpeciesProcedure(tmpSubstrate); // Update overall number of species and molecules
				reactionFlag = true;
            }else{
				ExitWithError("perform_endo_ComplexFormation", "Catalyst or Substrate are not avalaible");
            }
        }
        if(tmpNRGSide == CATALYSTLOAD)
        {
            if((allSpecies.at(tmpCatalyst).getChargeMols() > 0) && (allSpecies.at(tmpSubstrate).getNOTchargeMols() > 0))
            {
				unchargeMolProcess(tmpCatalyst);
				allSpecies.at(tmpCatalyst).decrement(volume);
				if(!allSpecies.at(tmpCatalyst).getConcentrationFixed())
					decMolSpeciesProcedure(tmpCatalyst); // Update overall number of species and molecules

				allSpecies.at(tmpSubstrate).decrement(volume);
				if(!allSpecies.at(tmpSubstrate).getConcentrationFixed())
					decMolSpeciesProcedure(tmpSubstrate); // Update overall number of species and molecules
				reactionFlag = true;
            }else{
				ExitWithError("perform_endo_ComplexFormation", "Catalyst or Substrate are not avalaible");
            }
        }
        if(tmpNRGSide == BOTHLOAD)
            {
            if((allSpecies.at(tmpCatalyst).getChargeMols() > 0) && (allSpecies.at(tmpSubstrate).getChargeMols() > 0))
            {
				unchargeMolProcess(tmpCatalyst);
				allSpecies.at(tmpCatalyst).decrement(volume);
				if(!allSpecies.at(tmpCatalyst).getConcentrationFixed())
					decMolSpeciesProcedure(tmpCatalyst); // Update overall number of species and molecules

				unchargeMolProcess(tmpSubstrate);
				allSpecies.at(tmpSubstrate).decrement(volume);
				if(!allSpecies.at(tmpSubstrate).getConcentrationFixed())
					decMolSpeciesProcedure(tmpSubstrate); // Update overall number of species and molecules
				reactionFlag = true;
            }else{
				ExitWithError("perform_endo_ComplexFormation", "Catalyst or Substrate are not avalaible");
            }
        }
	// EVALUATE COMPLEX (allSpecies.at(tmpCatalyst).getSequence().length() is definitly the cutting point)
	if(!complexEvaluation(tmpComplex, tmp__RndDoubleGen, allSpecies.at(tmpCatalyst).getSequence().length(),
                                                  tmpCatalyst, tmpSubstrate, tmpCatID, tmpSecSub, ESOERGONIC))
			ExitWithError("complexEvaluation", "Problems during the complex evalutation");
		
	incCpxFormCounter();

	if(debugLevel == FINDERRORDURINGRUNTIME) cout << "\t\tenvironment::perform_endo_ComplexFormation end" << endl;
	return reactionFlag;	
}


/**
 Perform COMPLEX DISASSOCIATION reaction
 @version 1.1
 @date 2010.06.08
 @param acs_longInt tmpComplex Complex ID
 @param acs_longInt tmpCatalyst Catalyst ID
 @param acs_longInt tmpSubstrate Substrate ID
 @param MTRand& tmp__RndDoubleGen random generator
 */
bool environment::performComplexDissociation(acs_longInt tmpComplex, acs_longInt tmpCatalyst, acs_longInt tmpSubstrate, MTRand& tmp__RndDoubleGen)
{
        if(debugLevel == FINDERRORDURINGRUNTIME) cout << "\t\tenvironment::performComplexDissociation start" << endl;

        if(debugLevel >= SMALL_DEBUG)
            cout << "\t\t\t|-COMPLEX " <<  tmpComplex << " will be degradated" << endl;
        bool reactionFlag = false;

    try{

        // Complex (loaded or not loaded) to dissociate is chosen according to the total amount of loaded and not loaded complexes present in the system
        acs_double ratioLoadedCpx = allSpecies.at(tmpComplex).getChargeMols() / allSpecies.at(tmpComplex).getAmount();
        acs_double rndNumb = tmp__RndDoubleGen();

        if(debugLevel >= SMALL_DEBUG)
        {
            cout << "\t\t\t\t|- BEFORE: Complexes[" << allSpecies.at(tmpComplex).getSequence() << "]- Loaded:" << allSpecies.at(tmpComplex).getChargeMols()
                 << "- Not Loaded:" << allSpecies.at(tmpComplex).getNOTchargeMols()
                 << " - Catalysts[" << allSpecies.at(tmpCatalyst).getSequence() << "]: " << allSpecies.at(tmpCatalyst).getAmount()
                 << " - Substrates[" << allSpecies.at(tmpSubstrate).getSequence() << "]: " << allSpecies.at(tmpSubstrate).getAmount() << endl;
        }

        if(rndNumb > ratioLoadedCpx)
        {
            // DECREMENT TOTAL AMOUNT OF THIS COMPLEX (not loaded)
            if(allSpecies.at(tmpComplex).getNOTchargeMols() < 1)
            {
                ExitWithError("performComplexDissociation", "There are not avalaible not loaded complexes to dissociate");
            }else{
                allSpecies.at(tmpComplex).decrement(volume);
                if(!allSpecies.at(tmpComplex).getConcentrationFixed())
                    decCpxProcedure(tmpComplex);
                reactionFlag = true;
            }
        }else{
            // DECREMENT TOTAL AMOUNT OF THIS COMPLEX (loaded)
            if(allSpecies.at(tmpComplex).getChargeMols() < 1)
            {
                ExitWithError("performComplexDissociation", "There are not avalaible loaded complexes to dissociate");
            }else{
                unchargeMolProcess(tmpComplex);
                allSpecies.at(tmpComplex).decrement(volume);
                if(!allSpecies.at(tmpComplex).getConcentrationFixed())
                    decCpxProcedure(tmpComplex);
                reactionFlag = true;
            }
        }

        // UPDATE CATALYST AND SUBSTRATE
        allSpecies.at(tmpCatalyst).increment(volume);
        if(!allSpecies.at(tmpCatalyst).getConcentrationFixed())
        {
            incMolProcedure(tmpCatalyst); // Increment total number of catalyst
            incSpeciesProcedure(tmpCatalyst);
        }

        allSpecies.at(tmpSubstrate).increment(volume);
        if(!allSpecies.at(tmpSubstrate).getConcentrationFixed())
        {
            incMolProcedure(tmpSubstrate); // Increment total number of substrate
            incSpeciesProcedure(tmpSubstrate);
        }

        if(debugLevel >= SMALL_DEBUG)
        {
            cout << "\t\t\t\t|- AFTER: Complexes[" << allSpecies.at(tmpComplex).getSequence() << "]- Loaded:" << allSpecies.at(tmpComplex).getChargeMols()
                 << "- Not Loaded:" << allSpecies.at(tmpComplex).getNOTchargeMols()
                 << " - Catalysts: " << allSpecies.at(tmpCatalyst).getAmount()
                 << " - Substrates: " << allSpecies.at(tmpSubstrate).getAmount() << endl;
        }

    }
    catch(exception&e)
    {
        cerr << "exceptioncaught:" << e.what() << endl;
        ExitWithError("something wrong in performComplexDissociation","exceptionerrorthrown");
    }

    incCpxDissCounter();

    if(debugLevel == FINDERRORDURINGRUNTIME) cout << "\t\tenvironment::performComplexDissociation end" << endl;

    return reactionFlag;
}

/**
 Perform SPONTANEOUS CLEAVAGE reaction
 @version 1.0
 @date 2013.10.28
 @param acs_longInt tmpReaction Reaction ID
 @param MTRand& tmp__RndDoubleGen random generator
 */
bool environment::performSpontaneousCleavage(acs_longInt tmpReaction, MTRand& tmp__RndDoubleGen)
{
        if(debugLevel == FINDERRORDURINGRUNTIME) cout << "\t\tenvironment::performComplexDissociation start" << endl;

        if(debugLevel >= SMALL_DEBUG)
            cout << "\t\t\t|-REACTION " <<  tmpReaction << " will spontaneously occur" << endl;
        bool reactionFlag = false;

        acs_longInt tmpSubstrate = allReactions.at(tmpReaction).getSpecies_I();
        acs_longInt tmpProduct_I = allReactions.at(tmpReaction).getSpecies_II();
        acs_longInt tmpProduct_II = allReactions.at(tmpReaction).getSpecies_III();

   try{

		// substrate decrement
		if(!allSpecies.at(tmpSubstrate).getAmount() > 0)
		{
			ExitWithError("performCleavage", "Substrate or catalyst not Avalaible!!!");
		}else{

			allSpecies.at(tmpSubstrate).decrement(volume);
			if(!allSpecies.at(tmpSubstrate).getConcentrationFixed())
				decMolSpeciesProcedure(tmpSubstrate); // Update overall number of species and molecules
			reactionFlag = true;

		}

	}catch(exception&e){

		cout<<"" << endl;
		cerr << "exceptioncaught:" << e.what() << endl;
		ExitWithError("something wrong in substrate decrement within performSpontaneousCleavage method","exceptionerrorthrown");
	}

	allReactions.at(tmpReaction).updateTotEvents();


	try{
		// EVALUATE PRODUCT 1
		if(!newSpeciesEvaluationIII(tmpProduct_I, tmp__RndDoubleGen))
			ExitWithError("newSpeciesEvaluationIII", "Problems during the new species evalutation");
	}catch(exception&e){
		cout<<" if(!newSpeciesEvaluation(allSpecies.at(tmpProduct_I).getSequence(), tmp__RndDoubleGen))" << endl;
		cerr << "exceptioncaught:" << e.what() << endl;
		ExitWithError("something wrong in evaluating product 1 in performCleavage method","exceptionerrorthrown");
	}

	try{
		// EVALUATE PRODUCT 2
		if(!newSpeciesEvaluationIII(tmpProduct_II, tmp__RndDoubleGen))
			ExitWithError("newSpeciesEvaluationIII", "Problems during the new species evalutation");
	}catch(exception&e){
		cout<<" if(!newSpeciesEvaluation(allSpecies.at(tmpProduct_II).getSequence(), tmp__RndDoubleGen))" << endl;
		cerr << "exceptioncaught:" << e.what() << endl;
		ExitWithError("something wrong in evaluating product 2 in performCleavage method","exceptionerrorthrown");
	}

	if(debugLevel == FINDERRORDURINGRUNTIME) cout << "\t\tenvironment::performCleavage end" << endl;
	// Update reaction counter
	incSpontDissCounter();
	return reactionFlag;
}

/**
 Perform SPONTANEOUS CONDENSATION reaction
 @version 1.0
 @date 2013-10-28
 @param acs_longInt tmpReaction Reaction ID
 @param MTRand& tmp__RndDoubleGen random generator
 */
bool environment::performSpontaneousCondensation(acs_longInt tmpReaction, MTRand& tmp__RndDoubleGen)
{
	if(debugLevel == FINDERRORDURINGRUNTIME) cout << "\t\tenvironment::performCondensation start" << endl;

	bool reactionFlag = false;
	acs_longInt tmpSub_I = allReactions.at(tmpReaction).getSpecies_II();
	acs_longInt tmpSub_II = allReactions.at(tmpReaction).getSpecies_III();
	acs_longInt tmpProduct = allReactions.at(tmpReaction).getSpecies_I();
	acs_longInt tmpSub_I_Amount;
	acs_longInt tmpSub_II_Amount;

    try{
		if(nrgBoolFlag == ENERGYBASED)
		{
			tmpSub_I_Amount = allSpecies.at(tmpSub_I).getNOTchargeMols();
			tmpSub_II_Amount = allSpecies.at(tmpSub_II).getNOTchargeMols();
		}else {
			tmpSub_I_Amount = allSpecies.at(tmpSub_I).getAmount();
			tmpSub_II_Amount = allSpecies.at(tmpSub_II).getAmount();
		}
    }catch(exception&e){
		cout<<"tmpSub_I_Amount = allSpecies.at(tmpSub_I).getNOTchargeMols(); " << endl;
		cerr << "exceptioncaught:" << e.what() << endl;
		ExitWithError("","exceptionerrorthrown");
    }


	// check whether all the species amounts to decrement are greater than 1
	if(tmpSub_I_Amount == 0)
		ExitWithError("performSpontaneousCondensation", "Substrate not Avalaible!!!");
	if(tmpSub_II_Amount == 0)
		ExitWithError("performSpontaneousCondensation", "Complex not Avalaible!!!");

	// The substrate and the complex are consumed (the first substrate has been already decremented during the complex formation reaction
	if ((tmpSub_I_Amount > 0) && (tmpSub_II_Amount > 0))
	{
		// REACTION!!!

        allSpecies.at(tmpSub_I).decrement(volume);
        if(!allSpecies.at(tmpSub_I).getConcentrationFixed())
            decMolSpeciesProcedure(tmpSub_I); // decrement total number of molecules and, if so, species

        allSpecies.at(tmpSub_II).decrement(volume);
        if(!allSpecies.at(tmpSub_II).getConcentrationFixed())
            decCpxProcedure(tmpSub_II); // decrement total number of complexes token and, if so, types

        reactionFlag = true;

	}else{
		ExitWithError("performSpontaneousCondensation", "Complex or second substrate are not avalaible");
	}
	// Update Reactions structure with a secure check on ID reaction
    allReactions.at(tmpReaction).updateTotEvents();

	// EVALUATE PRODUCT
	if(!newSpeciesEvaluationIII(tmpProduct, tmp__RndDoubleGen))
		ExitWithError("newSpeciesEvaluationIII", "Problems during the new species evalutation");
	if(debugLevel == FINDERRORDURINGRUNTIME) cout << "\t\tenvironment::performCondensation end" << endl;

	// Update reaction counter
	incSpontAssCounter();
	return reactionFlag;
}


//TR
/*
 Perform PHOSPHORILATION reaction
 @version 1.0
 @param acs_longInt tmpSpecies Species ID
 @date 2010-05-25

bool environment::performPhosphorilation(acs_longInt tmpSpecies)
{
	if(debugLevel == FINDERRORDURINGRUNTIME) cout << "\t\tenvironment::performPhosphorilation start" << endl;
	
	if(debugLevel >= SMALL_DEBUG)
		cout << "\t\t\t|-Species " << tmpSpecies << " will be Phosphorilated" << endl;
	bool reactionFlag = false;
	
	if(debugLevel >= SMALL_DEBUG)
	{
		cout << "\t\t\t\t|- BEFORE: loaded molecules[" << allSpecies[tmpSpecies].getSequence() << "]: " << allSpecies[tmpSpecies].getChargeMols() << endl;
		cout << "\t\t\t\t|- BEFORE: loaded carriers[" << nrgCarrier.getSequence() << "]: " << nrgCarrier.getChargeMols() << endl;
	}
	
	// Molecule phosphorilation
	if(!allSpecies.at(tmpSpecies).getNOTchargeMols() > 0)
	{
		ExitWithError("chargeMol", "Problems during the PHOSPHORILATION process, no molecules to load");
	}else{
		if(!allSpecies[tmpSpecies].chargeMol())
			ExitWithError("chargeMol", "Problems during the PHOSPHORILATION process, in particular loaded molecules are more than the unloaded molecules");
		else
			reactionFlag = true;
	}

        // Increase the counter of loaded molecules
        incOverallLoadedMolsCounter();

	
	// Decrease the number of loaded energy carriers
	if(!nrgCarrier.unchargeMol())
		ExitWithError("unchargeMol", "Problems during the PHOSPHORILATION process, in particular loaded carriers are more than the unloaded molecules");	
	else
		reactionFlag = true;
	
	if(debugLevel >= SMALL_DEBUG)
	{
		cout << "\t\t\t\t|- AFTER: loaded molecules[" << allSpecies[tmpSpecies].getSequence() << "]: " << allSpecies[tmpSpecies].getChargeMols() << endl;
		cout << "\t\t\t\t|- AFTER: loaded carriers[" << nrgCarrier.getSequence() << "]: " << nrgCarrier.getChargeMols() << endl;
	}
	
	if(debugLevel == FINDERRORDURINGRUNTIME) cout << "\t\tenvironment::performPhosphorilation end" << endl;
	
	return reactionFlag;
}*/

/**
 Evaluate new species (Optimized function with new species reactions update process)
 @version 3.0
 @param acs_int tmpNewSpecies New species ID to evaluate
 @param MTRand& tmp___RndDoubleGen random number generator
 */
bool environment::newSpeciesEvaluationIII(acs_longInt tmpNewSpecies, MTRand& tmp___RndDoubleGen)
{
        if(debugLevel == FINDERRORDURINGRUNTIME) cout << "\t\t\tenvironment::newSpeciesEvaluationIII start" << endl;

        bool newSpeciesFlag = true; // Function Flag Control
        bool tmpNotEqualSeqBetweenTwoSpecies = true; // True whether there is not equal species
        bool tmpAlreadyEvaluated = false; // True if the species has been already evaluated
        bool toEvaluate = true; // another control to check whether the evaluation is to do

        if(debugLevel == GILLESPIESTUFF) cout << "\t\t\t|- New Species: " << tmpNewSpecies << endl;

    try{
		if(allSpecies.at(tmpNewSpecies).getAmount() == 0)
			allSpecies.at(tmpNewSpecies).rebornsIncrement(); // IF THE SPECIES REAPPEAR THE REBORN COUNTER IS UPDATED
		allSpecies.at(tmpNewSpecies).increment(volume); // INCREMENT THE NUMBER OF ELEMENTS OF THIS SPECIES (conc fixed check is inside the function)
		if(!allSpecies.at(tmpNewSpecies).getConcentrationFixed())
		{
			incMolProcedure(tmpNewSpecies); // Increment overall amount of species
			incSpeciesProcedure(tmpNewSpecies);
		}

		tmpNotEqualSeqBetweenTwoSpecies = false; // IF THE SPECIES IS ALREADY PRESENT
		if(allSpecies.at(tmpNewSpecies).getEvaluated() == 1) // IF THE SPECIES HAS BEEN ALREADY EVALUATED TOO
		{
			tmpAlreadyEvaluated = true;
		}else{
			if(allSpecies.at(tmpNewSpecies).getComplexCutPnt() == 0) // IF THE SPECIES IS NOT A COMPLEX
			{
				// IF THE SPECIES IS ALREADY PRESENT BUT NOT EVALUATED
				// PRECIPITATION AND DIFFUSION CONTRIBUTE COMPUTATION
				acs_double tmpDiffusionContribute = createDiffusionRenforcement(diffusion_contribute, 1);
				allSpecies.at(tmpNewSpecies).setDiffusion(tmpDiffusionContribute);
				// IN THIS CASE SET SOLUBILITY REFERES TO THE SPECIES OBJECT
				bool tmpSolubilityState = setSolubility(allSpecies.at(tmpNewSpecies).getSequenceLength(), tmp___RndDoubleGen);
				allSpecies.at(tmpNewSpecies).setSolubility(tmpSolubilityState);
				allSpecies.at(tmpNewSpecies).setKphospho(K_nrg);
				allSpecies.at(tmpNewSpecies).setEvaluated(); // IF the reactions have not been yet created the molecules evaluation is setting on 1 and the
				// reactions will be created below
			}
		}


        //IF THE MOLECULES IS NOT YET EVALUATED
        if(tmpAlreadyEvaluated == false) //tmpAlreadyEvaluated = false at the beginning of the function
        {
            acs_longInt totalNumberOfConceivableReactions = 0;

            if(newSpeciesProbMinThreshold < ratioBetweenNewGillTotGill)
            {
				if(toEvaluate) // TO SILENT IF NEW SPECIES CREATION MUST BE AVOIDED!!!
				{
					// SET REACTIONS FOR THIS NEW SPECIES
					// COPY ALL THE ALREADY EVALUATED SPECIES ID IN A TEMPORARY VECTOR
					vector<acs_longInt> tmpAlreadyEvaluatedSpeciesVector;
					for(acs_longInt tmpS = 0; tmpS < (acs_longInt)allSpecies.size(); tmpS++)
					{
						if(allSpecies.at(tmpS).getEvaluated() == 1)
						{
							// DEEP AND SECURE CHECK
							if(allSpecies.at(tmpS).getComplexCutPnt() > 0)
							{
								ExitWithError("newSpeciesEvaluation", "A complex cannot be evaluated and involved in reactions!!!");
							}else{
								tmpAlreadyEvaluatedSpeciesVector.push_back(tmpS);
							}
						}
					}

					// COMPUTE TOTAL NUMBER OF CONCEIVABLE REACTIONS CONSIDERING ALREADY EVALUATED MOLECULES
					for(acs_int presSpeciID = 0; presSpeciID < (acs_longInt)tmpAlreadyEvaluatedSpeciesVector.size(); presSpeciID++)
					{
						totalNumberOfConceivableReactions += allSpecies.at(tmpAlreadyEvaluatedSpeciesVector.at(presSpeciID)).getSequenceLength() - 1;
					}
					totalNumberOfConceivableReactions += pow((double)tmpAlreadyEvaluatedSpeciesVector.size(), 2.0);

					if(debugLevel >= HIGH_DEBUG){printAllSpeciesIdAndSequence();}
					if(debugLevel >= RUNNING_VERSION)
					{
						cout << "\t\t|-------------------------------------------------------\\" << endl;
						cout << "\t\t|- !! EVENT !! New species has been evaluated" << endl;
						cout << "\t\t|- ID: " << allSpecies.at(tmpNewSpecies).getID() << " - Sequence: " << allSpecies.at(tmpNewSpecies).getSequence() << endl;
						cout << "\t\t|- Total new number of conceivable reactions:" << totalNumberOfConceivableReactions << endl;
					}

					// COMPUTE THE REAL NUMBER OF REACTIONS FOR THIS SPECIES
					acs_int reactionsForThisSpecies;
					reactionsForThisSpecies = computeSngSpeciesRctsNumber(totalNumberOfConceivableReactions, tmp___RndDoubleGen);

					if(debugLevel >= RUNNING_VERSION)
					{
						cout << "\t\t|- Number of reaction to catalyze: " << reactionsForThisSpecies << endl;
						cout << "\t\t|-------------------------------------------------------/" << endl;
					}

					// CREATE REACTIONS FOR THIS SPECIFIC SPECIES
					if((reactionsForThisSpecies > 0) && (allSpecies.at(tmpNewSpecies).getSolubility() == SOLUBLE))
					{
						if(!createReactionsForThisSpecies(tmpNewSpecies, reactionsForThisSpecies, tmp___RndDoubleGen, tmpAlreadyEvaluatedSpeciesVector, NEWREACTIONS))
							ExitWithError("createReactionsForThisSpecies", "	|- !*!*!*! Problem with the reactions creation");
					}

					if(debugLevel >= SMALL_DEBUG){cout << "\t\t|- DONE!!! " << endl;}

					// IF the new species is soluble
					if(allSpecies.at(tmpNewSpecies).getSolubility() == SOLUBLE)
					{
						// ------------------------------------------
						// OLD SPECIES REACTIONS UPDATE
						// ------------------------------------------
						// Compute theoretical number of reactions involving the new species.
						// CLEAVAGE reactions: allSpecies.at(tmpNewSpecies).getSequenceLength() - 1
						// CONDENSATION reactions: (2 * tmpAlreadyEvaluatedSpeciesVector.size())
						acs_int theoreticalNofNewRcts = (tmpAlreadyEvaluatedSpeciesVector.size() - 1) *
														((allSpecies.at(tmpNewSpecies).getSequenceLength() - 1) +
														(2 * tmpAlreadyEvaluatedSpeciesVector.size()));

						// According to new number of possibile new reactions, the probability and the number of species, the number of new reaction to catalyze is now computed
						// rctToAdd < theoreticalNofNewRcts
						acs_int rctToAdd = computeSngSpeciesRctsNumber(theoreticalNofNewRcts, tmp___RndDoubleGen);

						// For each new reaction
						acs_longInt speciesToUpdate;
						vector<acs_double> possRctType;
						possRctType.push_back(2 * tmpAlreadyEvaluatedSpeciesVector.size());
						possRctType.push_back(possRctType.at(0) + (allSpecies.at(tmpNewSpecies).getSequenceLength() - 1));
						acs_int rctType;

						if(debugLevel >= RUNNING_VERSION)
							cout << "\t\t\t|----- Old species catalysis have been update with " << rctToAdd << " new reactions ("
							<< theoreticalNofNewRcts << " conceivable reactions for all " << tmpAlreadyEvaluatedSpeciesVector.size() << " already eval species)" << endl;

						for(acs_int newR = 0; newR < rctToAdd; newR++)
						{
							// Select catalyst to update
							speciesToUpdate = tmpAlreadyEvaluatedSpeciesVector.at(getIntRandom(nonCatalyticLastID+1, tmpAlreadyEvaluatedSpeciesVector.size() - 1, tmp___RndDoubleGen));
							// According to the number of cleavage and condensation select wheater cleavage or condensation
							rctType = returnSelectionIdFromAWeightProbVector(possRctType, possRctType.at(1), tmp___RndDoubleGen, __LINE__);
							// According to the species to update and the reaction type new reactions are created
						    try{
						    	updateReactions(speciesToUpdate, tmpNewSpecies, rctType, tmpAlreadyEvaluatedSpeciesVector, tmp___RndDoubleGen);
						    }catch(exception&e)
						    {
						        cerr << "exceptioncaught:" << e.what() << endl;
						        ExitWithError("something wrong in update reactions method","exceptionerrorthrown");
						    }
						}
					} // if(allSpecies.at(tmpIdSpeciesToEvaluate).getSolubility() == SOLUBLE)
				} // end if to evaluate
            } // newSpeciesProbMinThreshold < ratioBetweenNewGillTotGill
        } // tmpAlreadyEvaluated == false

        if(debugLevel == FINDERRORDURINGRUNTIME) cout << "\t\t\tenvironment::newSpeciesEvaluationIII end" << endl;

    }
    catch(exception&e)
    {
        cerr << "exceptioncaught:" << e.what() << endl;
        ExitWithError("something wrong in newSpeciesEvaluationIII method","exceptionerrorthrown");
    }

    return newSpeciesFlag;
}
/**
 Complex evaluation
 @version 1.1
 @date 2010-06-04
 @param string tmpNewSpecies New species sequence to evaluate 
 @param MTRand& tmp___RndDoubleGen random number generator
 @param tmpCuttingPnt Complex cutting point
 @param tmpCatalyst_ID catalyst ID
 @param tmpSubstrate_ID substrate ID
 @param tmp_catalysisID catalysis ID
 @param tmpCpxType ENDOERGONIC or ESOERGONIC
 */
bool environment::complexEvaluation(string tmpComplex, MTRand& tmp___RndDoubleGen, acs_int tmpCuttingPnt, 
									acs_longInt tmpCatalyst_ID, acs_longInt tmpSubstrate_ID,
									acs_longInt tmpCatID, acs_longInt tmpSecSub_ID, bool tmpCpxType)
{
        if(debugLevel == FINDERRORDURINGRUNTIME) cout << "\t\t\tenvironment::ComplexEvaluation start" << endl;

        bool newComplexFlag = true; // Function Flag Control
        bool tmpNotEqualSeqBetweenTwoSpecies = true; // True whether there is not equal species
        acs_longInt tmpIdSpeciesToEvaluate = 0;

        if(debugLevel >= SMALL_DEBUG)
            cout << "\t\t\t|- New Species: " << tmpComplex << " CP: " << tmpCuttingPnt << endl;


    try{

        //CHECK WHETHER THE SPECIES (OR COMPLEX) IS ALREADY PRESENT
        for(acs_longInt i = 0; i < (acs_longInt)allSpecies.size(); i++)
        {
            if((allSpecies.at(i).getSequence() == tmpComplex) & (allSpecies.at(i).getComplexCutPnt() > 0)) // IF THE COMPLEX IS ALREADY PRESENT
            {
                if(allSpecies.at(i).getComplexCutPnt() == tmpCuttingPnt) // IF THE CUTTING POINT IS THE SAME
                {
                    allSpecies.at(i).increment(volume); // INCREMENT CHE NUMBER OF ELEMENTS OF THIS SPECIES THAT COULD BE ALSO A COMPLEX
                    if(allSpecies.at(i).getAmount() == 1){incNumberOfCpx();} // If the number of complex was 0 the total number of complex is increased
                    incNumberOfCpxMols(); // Increase the total number of complex copies

                    if(tmpCpxType == ESOERGONIC)
                    {
                        allSpecies.at(i).chargeMol();
                        incOverallLoadedMolsCounter();
                    } // If the complex is formed of loaded molecules it is considered loaded as well

                    tmpNotEqualSeqBetweenTwoSpecies = false; // IF THE SPECIES IS ALREADY PRESENT
                    tmpIdSpeciesToEvaluate = i; // SET THE ALREADY PRESENT ID OF THE SPECIES
                    //IF THE SPECIES IS A COMPLEX
                    // THE SPECIE IS A COMPLEX, ONLY DIFFUSION CONSTANT HAS TO BE CREATED, WE ASSUEME THAT DEGRADATION IS MORE FASTER THAN
                    // PRECIPITATION, HENCE PRECIPITATION IS SET TO 0 (NOT COMPUTED)
                    if(allSpecies.at(i).getDiffusionEnh() == 0)
                    {
                        acs_double tmpDiffusionContribute = createDiffusionRenforcement(diffusion_contribute, 1);
                        allSpecies.at(i).setDiffusion(tmpDiffusionContribute);
                    }
                    // Add the second substrate to the complex species
                    allSpecies.at(i).insertSecSub(tmpSecSub_ID,allCatalysis.at(tmpCatID).getKass(), tmpCatID);
                    break;
                }
            }
        }

        // IF THE COMPLEX IS ACTUALLY NEW
        if(tmpNotEqualSeqBetweenTwoSpecies == true)
        {
            // PRECIPITATION AND DIFFUSION CONTRIBUTE COMPUTATION
            acs_double tmpDiffusionContribute = createDiffusionRenforcement(diffusion_contribute, tmpComplex.length());
            // CREATE NEW SPECIES OBJECT AND STORE IT INTO THE ENVIRONMENT SPECIES VECTOR
            tmpIdSpeciesToEvaluate = (acs_longInt)allSpecies.size();


            // COMPLEX FORCED CREATIONS
            // ALSO tmp_catalysisID of this function IS TO REMOVE
            allSpecies.push_back(species(tmpIdSpeciesToEvaluate, tmpComplex, tmpDiffusionContribute,
                                         SOLUBLE, K_cpxDiss, tmpCuttingPnt,
                                         tmp___RndDoubleGen, tmpCatalyst_ID, tmpSubstrate_ID, volume, 0, NOTENERGIZABLE));
            incNumberOfCpx();	// Increment number of complex types
            incNumberOfCpxMols();	// Increment number of complex tokens

            //if(tmpIdSpeciesToEvaluate == 247)
            	//printAllSpeciesIdAndSequence();

            // Add the second substrate to the complex species
            if(debugLevel == COMPLEXSTUFF) cout << "COMPLEX STUFF |- Catalysis: " << tmpCatID << endl;
            allSpecies.at(allSpecies.size()-1).insertSecSub(tmpSecSub_ID,allCatalysis.at(tmpCatID).getKass(), tmpCatID);

            if(tmpCpxType == ESOERGONIC)
            {
                allSpecies.at(tmpIdSpeciesToEvaluate).chargeMol();
                incOverallLoadedMolsCounter();
            } // If the complex is formed by a loaded molecules, it is loaded as well

            if(debugLevel >= SMALL_DEBUG)
                cout << "\t\t\t\t|- New Complex " <<  tmpComplex << " has been created" << endl;
        }

    }
    catch(exception&e)
    {
        cerr << "exceptioncaught:" << e.what() << endl;
        ExitWithError("error in complexEvaluation method","exceptionerrorthrown");
    }

    if(debugLevel == FINDERRORDURINGRUNTIME) cout << "\t\t\tenvironment::ComplexEvaluation end" << endl;

    return newComplexFlag;
}

/**
 Update the species age
 @version 1.0
 */
void environment::updateSpeciesAges()
{
    try{
        if((acs_longInt)allSpecies.size() > 0)
        {
            for(vector<species>::iterator tmpAllSpeciesIter = allSpecies.begin(); tmpAllSpeciesIter != allSpecies.end(); tmpAllSpeciesIter++)
            {
                // UPDATE SPECIES AGE ----------------
                if((tmpAllSpeciesIter->getAmount() > 0) & (tmpAllSpeciesIter->getComplexCutPnt() == 0))
                {
                    tmpAllSpeciesIter->setNewAge(reactionsTime.at((reactionsTime.size() - 1)) - reactionsTime.at((reactionsTime.size() - 2)));
                }else{
                    tmpAllSpeciesIter->setNewAge(0); // IF THE SPECIES IS DEAD THE AGE OF LAST DEATH REMAINS
                }

            }
        }

    }
    catch(exception&e)
    {
        cerr << "exceptioncaught:" << e.what() << endl;
        ExitWithError("error in updateSpeciesAges method","exceptionerrorthrown");
    }
}

/**
 Change volume function
 @versione 1.0
 @date 2013/07/17
 @author Alessandro Filisetti
 */
void environment::changeVolume(acs_int tmpTimeSinceLastReaction)
{
	//volume = volume + ((volume * 10) * tmpTimeSinceLastReaction) ;
	cout << "volue change" << endl;
	volume = volume * 10 ;
}

/* ************************
 |
 |	PROMPT FUNCTIONS
 |
 *************************/

/**
 Shows all parameters uploaded from the configuration file
 @version 1.0
 */
void environment::showGlobalParameter()
{
	if(debugLevel == FINDERRORDURINGRUNTIME) cout << "environment::showGlobalParameter start" << endl;

	if(debugLevel >= RUNNING_VERSION)
	{ 
        cout << endl << endl;
        cout << "|   ****         ***   *     *          *****   *****                             " << endl;
        cout << "|  *            *   *  **    *         *       *                                  " << endl;
        cout << "| *        ***  *   *  * *   *  ****   *       *                                  " << endl;
        cout << "| *          *  ****   *  *  *  *  *    ****    ****       " << endl;
        cout << "| *       ****  **     *   * *  ****        *       *                             " << endl;
        cout << "|  *      *  *  * *    *    **  *           *       *                             "  << endl;
        cout << "|   ****  ****  *  *   *     *   ***   *****   *****                              " << endl;
        cout << "|  												                               " << endl;
        cout << "| CAtalytic     Reaction NEtworks      Stochastic Simulator	                   " << endl;
        cout << "| VERSION_ " << __SOFTVERSION__ << endl;
        cout << "| Author_ Alessandro Filisetti" << endl;
        cout << "| mail_ alessandro.filisetti@gmail.com" << endl;
        cout << endl << endl << endl;

		//TR if(onlyEnvironmentCreation)
		//TR 	out << "\t|-ONLY THE ENVIRONMENT WILL BE CREATED..." << endl; 
        cout << "|- GLOBAL PARAMETERS" << endl;
        cout << "\t|- Number of Generations: " << nGEN << endl;
        cout << "\t|- Number of Simulations: " << nSIM << endl;
        cout << "\t|- Number of Seconds: " << (double)nSeconds << endl;
        cout << "\t|- Number of Reactions Permitted: " << nReactions << endl;
        cout << "\t|- Max number of attempts in simulating the same network with different random seeds: " << nAttempts << endl;
        cout << "\t|- Max number of hours of the simulation (computational time): " << (double)nHours << endl;
        cout << "\t|- Last Firing disk species ID (in structure files upload configuration): " << lastFiringDiskSpeciesID << endl;
        cout << "\t|- time Structures Saving Interval: " << (double)timeStructuresSavingInterval << endl;
        cout << "\t|- Minimal new species prob to allow system expansion: " << (double)newSpeciesProbMinThreshold << endl;
        cout << "\t|- file times Saving Interval: " << (double)fileTimesSaveInterval << endl;
        cout << "\t|- Max lenght of non catalytic species: " << nonCatalyticMaxLength << endl;
		cout << "\t|- Overall initial concentration: " << overallConcentration << endl;
		cout << "\t|- Energy Carriers concentration: " << ECConcentration << endl;
		cout << "\t|- Reaction Probability: " << reactionProbability << endl;
		cout << "\t|- cleavage Probability: " << cleavageProbability << endl;
		cout << "\t|- Reverse reactions as main reactions possibility: " << main_rev_rct_allowed << endl;
		cout << "\t|- Reverse Reactions: " << reverseReactions << endl;
        cout << "\t|- Ratio between forward and backward reactions: " << revRctRatio << endl;
		cout << "\t|- ASSOCIATION kinetic constant: " << K_ass << endl;
		cout << "\t|- DISSOCIATION kinetic constant: " << K_diss << endl;
		cout << "\t|- COMPLEX FORMATION kinetic constant: " << K_cpx << endl;
		cout << "\t|- COMPLEX dissociation : " << K_cpxDiss << endl;
		cout << "\t|- ENERGY PHOSFORILATION constant: " << K_nrg << endl;
        cout << "\t|- ENERGY DECAY CONSTANT: " << K_nrg_decay << endl;
        cout << "\t|- SPONTANOUES CONDENSATION kinetic constant: " << K_spont_ass << endl;
        cout << "\t|- SPONTANOUES CLEAVAGE kinetic constant: " << K_spont_diss << endl;
		cout << "\t|- SOLUBILITY THRESHOLD: " << solubility_threshold << endl;
        if(influx_rate > 0)
		{
            cout << "\t|- The system is OPEN" << endl;
			cout << "\t|- Max out length: " << maxLOut << endl;
			cout << "\t|- MOLECULE DECAY kinetic constant: " << moleculeDecay_KineticConstant << endl;
			cout << "\t|- Refill rate: " << influx_rate << " conc/sec" << endl;
		}else{
			cout << "\t|- The system is CLOSED" << endl;
		}
        if(spontRct > 0){cout << "\t|- Spontaneuos Reactions are present" << endl;}
        else{cout << "\t|- Spontanoues Reactions are NOT present" << endl;}

		cout << "\t|- Alphabet: " << alphabet << endl;
		cout << "\t|- Volume: "<< volume << endl;
		cout << "\t|- Volume Growth: " << volumeGrowth << endl;
		cout << "\t|- Stochastic Division: " << stochDivision << endl;
		cout << "\t|- Random Seed: " << randomSeed << endl;	
		cout << "\t|- NEPERO: " << NEP << endl;
		cout << "\t|- AVOGADRO: " << AVO << endl;
		cout << "\t|- ENERGY: " << energy << endl;
		cout << "\t|- ratioSpeciesEnergizable: " << ratioSpeciesEnergizable << endl << endl;
		//TR cout << "\t|- Percentage of laoded carriers in the influx: " << percLoadedECInflux << endl;
		//TR switch (energyTarget) {
		//	case SUBSTRATELOAD:
		//		cout << "\t|- energyTarget: SUBSTRATES have to be loaded in endoergonic reactions" << endl;
		//		break;
		//	case CATALYSTLOAD:
		//		cout << "\t|- energyTarget: CATALYSTS have to be loaded in endoergonic reactions" << endl;
		//		break;
		//	default:
		//		cout << "\t|- energyTarget: BOTH SUBSTRATES AND CATALYSTS have to be loaded in endoergonic reactions" << endl;
		//	break;
		//}
	}	
	if(debugLevel == FINDERRORDURINGRUNTIME) cout << "environment::showGlobalParameter end" << endl;

}

/**
 Show all initial species in table format
 @version 1.0
*/
void environment::printInitialCondition()
{
	if(debugLevel == FINDERRORDURINGRUNTIME) cout << "environment::printInitialCondition start" << endl;

		cout << "\tNumber of theoretical species after the creation process -> " << this->getNumberOfTheoreticalSpecies() << endl;
		
		cout << "" << endl;
		cout << "ID" << "\t" << "SEQ" << "\t" << "N" << "\t" << "DIFF" << "\t" << "PREC" << "\t" << "DEG" << "\t" << "CUT" <<endl;
		
		for(acs_longInt i = 0; i < getNumberOfTheoreticalSpecies(); i++)
		{
			cout << allSpecies.at(i).getID() << "\t" << allSpecies.at(i).getSequence() << "\t" << 
				allSpecies.at(i).getAmount() << "\t" << allSpecies.at(i).getDiffusionEnh() << "\t" << 
				allSpecies.at(i).getSolubility() << "\t" << allSpecies.at(i).getComplexDegEnh() << "\t" <<
				allSpecies.at(i).getComplexCutPnt() << endl;
		}
		
		cout << "" << endl;	
		cout << " |-----------|" << endl;
		cout << " | REACTIONS |" << endl;
		cout << " |-----------|" << endl;	
		cout << "" << endl;
		if(getNumberOfReactions() > 0)
		{
			for(acs_longInt i = 0; i < getNumberOfReactions(); i++)
			{
				cout << allReactions.at(i).getID() << "\t" << allReactions.at(i).getType() << "\t" << 
						allReactions.at(i).getSpecies_I() << "\t" << allReactions.at(i).getSpecies_II() << "\t" << 
						allReactions.at(i).getSpecies_III() << endl;
			}
		}else{
			ExitWithError("printInitialCondition", "!!! THERE ARE NO REACTIONS IN THIS ENVIRONMENT");
		}
	
	cout << "" << endl;	
	cout << " |-----------|" << endl;
	cout << " | CATALYSIS |" << endl;
	cout << " |-----------|" << endl;	
	cout << "" << endl;
	if(getNumberOfCatalysis() > 0)
	{
		for(acs_longInt i = 0; i < getNumberOfCatalysis(); i++)
		{
			cout << allCatalysis.at(i).getCatId() << "\t" << allCatalysis.at(i).getCat() << "\t" << 
					allCatalysis.at(i).getReactionID() << "\t" << allCatalysis.at(i).getTotAmount()
					 << "\t" << allCatalysis.at(i).getCpxTarget() << endl;
		}
	}else{
		ExitWithError("printInitialCondition", "!!! THERE ARE NO CATALYSED REACTIONS!!!");
	}
	if(debugLevel == FINDERRORDURINGRUNTIME) cout << "environment::printInitialCondition end" << endl;

}

/**
 Show all the species with their ID
 @version 1.0
 */
void environment::printAllSpeciesIdAndSequence()
{
	if(debugLevel == FINDERRORDURINGRUNTIME) cout << "environment::printAllSpeciesIdAndSequence start" << endl;

	cout << "ID\tS" << endl;
	for(acs_longInt i = 0; i < (acs_longInt)allSpecies.size(); i++)
	{
		cout << allSpecies.at(i).getID() << "\t" << allSpecies.at(i).getSequence() << "\t" << allSpecies.at(i).getConcentration() 
			 << "\t" << allSpecies.at(i).getAmount() 
			 << "\t" << allSpecies.at(i).getSequence()
			 << "\t" << allSpecies.at(i).getComplexDegEnh()
			 << "\t" << allSpecies.at(i).getComplexCutPnt()
			 << "\t" << allSpecies.at(i).getEvaluated()
			 << "\t" << allSpecies.at(i).getCatalyst_ID()
			 << "\t" << allSpecies.at(i).getSubstrate_ID()
			 //<< "\t" << allSpecies.at(i).getNOTchargeMols()
			 //<< "\t" << allSpecies.at(i).getChargeMols()

			 << endl;
	}
	if(debugLevel == FINDERRORDURINGRUNTIME) cout << "environment::printAllSpeciesIdAndSequence end" << endl;

}

/**
 Show all the Gillespie Structure
 @version 1.0
 */
void environment::printGillespieStructure()
{
	if(debugLevel == FINDERRORDURINGRUNTIME) cout << "environment::printGillespieStructure start" << endl;

	cout << "---- GILLESPIE STRUCTURE ---- (SIZE " << allGillespieScores.size() << ")"<< endl << endl;
	cout << "ID\tSCORE\t\tRCT\tM1\tM2\tM3\tM4\tRCT\tCAT" << endl;
	for(acs_longInt i = 0; i < (acs_longInt)allGillespieScores.size(); i++)
	{
		cout << allGillespieScores.at(i).getID() << "\t" <<
			    allGillespieScores.at(i).getScore() << "\t\t" <<
				allGillespieScores.at(i).getIdReactionType() << "\t" <<
				allGillespieScores.at(i).getMolI() << "\t" <<
				allGillespieScores.at(i).getMolII() << "\t" <<
				allGillespieScores.at(i).getMolIII() << "\t" <<
				allGillespieScores.at(i).getMolIV() << "\t" <<
				allGillespieScores.at(i).getReactionID() << "\t" <<
				allGillespieScores.at(i).getCatalysisID() << endl;
	}
	if(debugLevel == FINDERRORDURINGRUNTIME) cout << "environment::printGillespieStructure end" << endl;

}

/**
 Print all the energized molecules
 @versione 1.0
 @date 2010-11-14
 */

void environment::printAllChargeMols()
{
	acs_int counter = 0;
	cout << "-------------" << endl;
        cout << "ID" << "\t" << "TOT" << "\t" << "-" << "\t" << "cum" << endl;
    for(vector<acs_double>::iterator tmpAllUnchargedIter = cumUnchargedAmountList.begin(); tmpAllUnchargedIter != cumUnchargedAmountList.end(); tmpAllUnchargedIter++)
	{
		cout << unchargedIDlist.at(counter) << "\t" << allSpecies.at(unchargedIDlist.at(counter)).getAmount() << "\t" 
			 << allSpecies.at(unchargedIDlist.at(counter)).getNOTchargeMols() << "\t" << *tmpAllUnchargedIter << endl;

                if(counter > 0)
                {
                    if((allSpecies.at(unchargedIDlist.at(counter)).getNOTchargeMols()) !=
                        (cumUnchargedAmountList.at(counter) - cumUnchargedAmountList.at(counter-1)))
                             ExitWithError("printAllChargeMols", "	|- !*!*!*! Cumulative sum is not right");
                }

		counter++;
	}
}
/**
 Clear all structures after each simulation
 @version 1.0
 */
void environment::clearAllStructures()
{
	if(debugLevel >= SMALL_DEBUG)
	{
		cout << "\t\t|- Number of Species: " << allSpecies.size() << endl;
		cout << "\t\t|- Number of Reactions: " << allReactions.size() << endl;
		cout << "\t\t|- Number of Catalysis: " << allCatalysis.size() << endl;
		cout << "\t\t|- Number of Times: " << allTimes.size() << endl;
		cout << "\t\t|- Number of Reactions Time: " << reactionsTime.size() << endl;
		cout << "\t\t|- Number of Reactions Selected: " << gillespieReactionsSelected.size() << endl;		
	}
       // printAllSpeciesIdAndSequence();
	allSpecies.clear();
	allReactions.clear();
	allCatalysis.clear();
	allTimes.clear();
	reactionsTime.clear();
	gillespieReactionsSelected.clear();
	if(debugLevel >= SMALL_DEBUG)
	{
		cout << "\t\t|- Number of Species: " << allSpecies.size() << endl;
		cout << "\t\t|- Number of Reactions: " << allReactions.size() << endl;
		cout << "\t\t|- Number of Catalysis: " << allCatalysis.size() << endl;
		cout << "\t\t|- Number of Times: " << allTimes.size() << endl;
		cout << "\t\t|- Number of Reactions Time: " << reactionsTime.size() << endl;
		cout << "\t\t|- Number of Reactions Selected: " << gillespieReactionsSelected.size() << endl;		
	}
	allSpecies = initialAllSpecies;
        //printAllSpeciesIdAndSequence();
	allReactions = initiaAllReactions;
	allCatalysis = initialAllCatalysis;
	

	if(debugLevel >= SMALL_DEBUG)
	{
		cout << "\t\t|- Number of Species: " << allSpecies.size() << endl;
		cout << "\t\t|- Number of Reactions: " << allReactions.size() << endl;
		cout << "\t\t|- Number of Catalysis: " << allCatalysis.size() << endl;
		cout << "\t\t|- Number of Times: " << allTimes.size() << endl;
		cout << "\t\t|- Number of Reactions Time: " << reactionsTime.size() << endl;
		cout << "\t\t|- Number of Reactions Selected: " << gillespieReactionsSelected.size() << endl;		
	}
	
	// Reset total amount of species, complexes and energy carriers for the next simulation

	numberOfSpecies = getTotalNumberOfSpecies(); // Set the initial number of species equal to the firing disk size
	numberOfNewSpecies = 0;
	numberOfMolecules = getTotalNumberOfMolecules();
	numberOfNewMolecules = 0;
	numberOfCpx = getTotalNumberOfComplexSpecies();
	numberOfCpxMols = getTotalNumberOfComplexes();
	timeSinceTheLastInFlux = 0;
	decimalMoleculesToEfflux = 0;
	decimalMoleculesToLoad = 0;
	decimalMoleculesToUNLOAD = 0;
	decimalComplexesToDissociate = 0;
    internalTimesStoredCounter = 0;
	setActualTime(0);
	resetReactionsCounter();
	lastEvaluatedSpeceisForNewReactions = 0;
	
}

/**
 set the concentrations to the initial values and reset internal statistics and counter
 @version 1.0
 */
void environment::resetConcentrationToInitialConditions()
{
    if(debugLevel >= RUNNING_VERSION)
    {
            cout << endl << "\t\t|- OLD STRUCTURES" << endl;
            cout << "\t\t|- Number of Species: " << allSpecies.size() << endl;
            cout << "\t\t|- Number of Reactions: " << allReactions.size() << endl;
            cout << "\t\t|- Number of Catalysis: " << allCatalysis.size() << endl;
            cout << "\t\t|- Number of Times: " << allTimes.size() << endl;
            cout << "\t\t|- Number of Reactions Time: " << reactionsTime.size() << endl;
            cout << "\t\t|- Number of Reactions Selected: " << gillespieReactionsSelected.size() << endl;
    }

    // RESET REACTIONS COUNTER
    for(vector<reactions>::iterator tmpAllReactionsIter = allReactions.begin(); tmpAllReactionsIter != allReactions.end(); tmpAllReactionsIter++)
    {
        tmpAllReactionsIter->resetEventsCounter();
    }
    // RESET CATALYSIS COUNTER
    for(vector<catalysis>::iterator tmpAllCatalysisIter = allCatalysis.begin(); tmpAllCatalysisIter != allCatalysis.end(); tmpAllCatalysisIter++)
    {
        tmpAllCatalysisIter->resetEventsCounter();
    }
    // RESET SPECIES PROPRIETIES
    for(vector<species>::iterator tmpAllSpeciesIter = allSpecies.begin(); tmpAllSpeciesIter != allSpecies.end(); tmpAllSpeciesIter++)
    {
        tmpAllSpeciesIter->resetAge();
        tmpAllSpeciesIter->resetReborns();
        tmpAllSpeciesIter->resetToInitConc(volume);
    }

    //allSpecies.clear();
    //allReactions.clear();
    //allCatalysis.clear();
	allTimes.clear();
	reactionsTime.clear();
	gillespieReactionsSelected.clear();

	if(debugLevel >= RUNNING_VERSION)
        {
            cout << "\t\t|- Number of Times: " << allTimes.size() << endl;
            cout << "\t\t|- Number of Reactions Time: " << reactionsTime.size() << endl;
            cout << "\t\t|- Number of Reactions Selected: " << gillespieReactionsSelected.size() << endl;
	}

	//LOADING SPECIES FROM FILE
//	if(!createInitialMoleculesPopulationFromSpecificFile(tmpSpeciesFilePath, tmpActGEN, tmpActSIM))
//		ExitWithError("createInitialMoleculesPopulationFromSpecificFile", "Problem with the species loading process");
//	//LOADING REACTIONS AND CATALYSIS FROM FILE
//	if(!createInitialReactionsLayerFromSpecificFile(tmpSpeciesFilePath, tmpActGEN, tmpActSIM))
//		ExitWithError("createInitialReactionsLayerFromSpecificFile", "Problem with the reactions loading process");
//	if(!createInitialCatalysisLayerFromSpecificFile(tmpSpeciesFilePath, tmpActGEN, tmpActSIM))
//		ExitWithError("createInitialCatalysisLayerFromSpecificFile", "Problem with the catalysis loading process");
	//CHECK STRUCTURES
	if(!structureCoherenceCheckUp())
	{
		cout << endl;
		ExitWithError("structureCoherenceCheckUp", "PROBLEM WITH STRUCTURE COHERENCE... BE CARE TO THE SIMULATION OUTCOMES!!!");
	}
	//STORE INITIAL STRUCTURES DATA IN ORDER TO REINITIALIZE THE STRUCTURE AFTER EACH SIMULATION
	storeInitialStructures();
	
	if(debugLevel >= RUNNING_VERSION)
	{
		cout << endl << "\t\t|- NEW STRUCTURES" << endl;
		cout << "\t\t|- Number of Species: " << allSpecies.size() << endl;
		cout << "\t\t|- Number of Reactions: " << allReactions.size() << endl;
		cout << "\t\t|- Number of Catalysis: " << allCatalysis.size() << endl;
		cout << "\t\t|- Number of Times: " << allTimes.size() << endl;
		cout << "\t\t|- Number of Reactions Time: " << reactionsTime.size() << endl;
		cout << "\t\t|- Number of Reactions Selected: " << gillespieReactionsSelected.size() << endl;		
	}

        // Reset total amount of species, complexes and energy carriers for the next simulation

	numberOfSpecies = getTotalNumberOfSpecies(); // Set the initial number of species equal to the firing disk size
	numberOfNewSpecies = 0;
	numberOfMolecules = getTotalNumberOfMolecules();
	numberOfNewMolecules = 0;
	numberOfCpx = getTotalNumberOfComplexSpecies();
	numberOfCpxMols = getTotalNumberOfComplexes();
	timeSinceTheLastInFlux = 0;
	decimalMoleculesToEfflux = 0;
	decimalMoleculesToLoad = 0;
	decimalMoleculesToUNLOAD = 0;
	decimalComplexesToDissociate = 0;
	setActualTime(0);
	resetReactionsCounter();
	
}

/**
 Store initial structures into storing variables
 @version 1.0
 */
void environment::storeInitialStructures()
{
	initialAllSpecies.clear();
	initiaAllReactions.clear();
	initialAllCatalysis.clear();
	initialAllSpecies = allSpecies;
	initiaAllReactions = allReactions;
	initialAllCatalysis = allCatalysis;
}

/* ************************
 |
 |	SAVE TO FILE FUNCTIONS
 |
 *************************/

/**
 This function creates a chain of zero as STRING according to tmpTotN and tmpCurrent N in order to make possible a sorting
 (e.g. tmpTotN = 1000, tmpCurrentN = 3, return 0003
 @version 1.0
 @param acs_int tmpTotN Total N
 @param acs_int tmpCurrentN current N
 */
string environment::zeroBeforeStringNumberSTD(acs_int tmpTotN, acs_int tmpCurrentN)
{
    string strZeroReturned = "";
    stringstream tempStrTmpTotN;
    tempStrTmpTotN << tmpTotN;
    stringstream tempStrtmpCurrentN;
    tempStrtmpCurrentN << tmpCurrentN;
    string QStrTmpTotN = tempStrTmpTotN.str();
    string QStrtmpCurrentN = tempStrtmpCurrentN.str();

    for(int i = 0; (unsigned)i < (unsigned)QStrTmpTotN.size() - (unsigned)QStrtmpCurrentN.size(); i++)
    {
        strZeroReturned += "0";
    }

    return strZeroReturned;
}

/**
 Save a file with the configuration parameters
 @version 1.0
 @date 2013/07/03
 */
bool environment::saveConfigurationFileSTD(string tmpStoringPath)
{
    if(debugLevel == FINDERRORDURINGRUNTIME) cout << "environment::saveConfigurationFileSTD start" << endl;

    cout << "Saving configuration file... ";
    string strConfigurationFile = tmpStoringPath + "/acsm2s.conf";

    ofstream fidFile;

    try{
        fidFile.open(strConfigurationFile.c_str(), ios::out | ios::app);
    }
    catch(exception&e)
    {
        cerr << "exceptioncaught:" << e.what() << endl;
        ExitWithError("error in method saveInfluxStructureSTD","exceptionerrorthrown");
    }

    fidFile << "# =================" << endl;
    fidFile << "# ACSM2S PARAMETERS" << endl;
    fidFile << "# Ver. " << __SOFTVERSION__ << endl;
    fidFile << "# =================" << endl << endl;

    fidFile << "# -----------------" << endl;
    fidFile << "# SYSTEM PARAMETERS" << endl;
    fidFile << "# -----------------" << endl << endl;

    //TR out << "# Limit the work of the software to the environment creation" << endl;
    //TR out << "onlyEnvironmentCreation=" << onlyEnvironmentCreation << endl << endl;

    fidFile << "# Number of Generations" << endl;
    fidFile << "nGEN=" << nGEN << endl << endl;

    fidFile << "# Number of Simulations per generation (after noMultipleSimsGens generations)" << endl;
    fidFile << "nSIM=" << nSIM << endl << endl;

    fidFile << "# Number of seconds" << endl;
    fidFile << "nSeconds=" << (double)nSeconds << endl << endl;

    fidFile << "# Number of max reactions permitted (step)" << endl;
    fidFile << "nReactions=" << nReactions << endl << endl;

    fidFile << "# Max number of hours (computational time) of the simulation (if 0 no limits are set)" << endl;
    fidFile << "nHours=" << (double)nHours << endl << endl;

    fidFile << "# Number of attempts in same network / different random seed (if 0 no limits are set) " << endl;
    fidFile << "nAttempts=" << nAttempts << endl << endl;

    fidFile << "# random seed (random if 0)" << endl;
    fidFile << "randomSeed=" << (int)randomSeed << endl << endl;

    fidFile << "# Debug Level: Runtime number of prompted messages level" << endl;
    fidFile << "debugLevel=" << debugLevel << endl << endl;

    fidFile << "# Save structures to file every..." << endl;
    fidFile << "timeStructuresSavingInterval=" << (double)timeStructuresSavingInterval << endl << endl;

    fidFile << "# Save file times avery..." << endl;
    fidFile << "fileTimesSaveInterval=" << (double)fileTimesSaveInterval << endl << endl;

    fidFile << "# Minimal new species creation probability to allow system expansion" << endl;
    fidFile << "# (If you want to avoid the creation of both new species and reactions set >= 1)" << endl;
    fidFile << "newSpeciesProbMinThreshold=" << (double)newSpeciesProbMinThreshold << endl << endl;

    fidFile << "# ------------------------" << endl;
    fidFile << "# ENVIRONMENTAL PARAMETERS" << endl;
    fidFile << "# ------------------------" << endl << endl;

    fidFile << "# Total number of species in the initial population" << endl;
    fidFile << "# 	If The number is greater than the theoretical number" << endl;
    fidFile << "#	according to initial max length and the alphabet all species" << endl;
    fidFile << "#	will be created, otherwise if it is equal to 0 the system will be" << endl;
    fidFile << "#	upload from file." << endl << endl;

    fidFile << "# Identificator of the last firing disk species" << endl;
    fidFile << "lastFiringDiskSpeciesID=" << lastFiringDiskSpeciesID << endl << endl;

    fidFile << "# Initial distribution" << endl;
    fidFile << "# 1- Proportional: Same number of molecules for each species" << endl;
    fidFile << "# 2- Uniform: Same number of molecules for each length" << endl;
    fidFile << "# 3- inversely proportional: more morecules for the short species" << endl << endl;

    fidFile << "# Overall initial concentration" << endl;
    fidFile << "overallConcentration=" << (double)overallConcentration << endl << endl;

    fidFile << "# Overall energy carriers concentration" << endl;
    fidFile << "ECConcentration=" << (double)ECConcentration << endl << endl;

    fidFile << "# Alphabet (e.g. AGCT for DNA)" << endl;
    fidFile << "alphabet=" << alphabet.c_str() << endl << endl;

    fidFile << "# Volume (dm^3)" << endl;
    fidFile << "volume=" << (double)volume << endl << endl;

    fidFile << "# Volume growth possibility (1: volume changes, 0: volume is fixed)" << endl;
    fidFile << "volumeGrowth=" << volumeGrowth << endl << endl;

    fidFile << "# Division type (1: Stochastic, 0: deterministic)" << endl;
    fidFile << "stochDivision=" << stochDivision << endl << endl;

    fidFile << "# ------------------" << endl;
    fidFile << "# DYNAMIC PARAMETERS" << endl;
    fidFile << "# ------------------" << endl << endl;

    fidFile << "# Energy introduction" << endl;
    fidFile << "#	0: Not energy in the system" << endl;
    fidFile << "#	1: Energy is present" << endl;
    fidFile << "energy=" << (double)energy << endl << endl;

    fidFile << "# Ratio of energizable species" << endl;
    fidFile << "ratioSpeciesEnergizable=" << (double)ratioSpeciesEnergizable << endl << endl;

    fidFile << "# IF 1 also monomers can catalyze reactions, otherwise reactions are catalyzed" << endl;
    fidFile << "# starting from dimers" << endl;
    fidFile << "nonCatalyticMaxLength=" << nonCatalyticMaxLength << endl << endl;

    fidFile << "# Catalysis probability" << endl;
    fidFile << "reactionProbability=" << (double)reactionProbability << endl << endl;

    fidFile << "# Cleavage probability (Condensation Probability is 1 - cleavage probability)" << endl;
    fidFile << "cleavageProbability=" << (double)cleavageProbability << endl << endl;

    fidFile << "# Enable the possibility of having reactions AS MAIN REACTIONS" << endl;
    fidFile << "main_rev_rct_allowed=" << main_rev_rct_allowed << endl << endl;

    fidFile << "# Enable reverse reactions" << endl;
    fidFile << "reverseReactions=" << reverseReactions << endl << endl;

    fidFile << "# Ratio between forward and backward reactions (if reverseReactions = TRUE)" << endl;
    fidFile << "revRctRatio=" << (double)revRctRatio << endl << endl;

    fidFile << "# Spontaneous Reactions (if present = TRUE)" << endl;
    fidFile << "spontRct=" << (double)spontRct << endl << endl;

    fidFile << "# kinetic constants" << endl;
    fidFile << "K_ass=" << (double)K_ass << endl;
    fidFile << "K_diss=" << (double)K_diss << endl;
    fidFile << "K_cpx=" << (double)K_cpx << endl;
    fidFile << "K_cpxDiss=" << (double)K_cpxDiss << endl;
    fidFile << "K_nrg=" << (double)K_nrg << endl;
    fidFile << "K_nrg_decay=" << (double)K_nrg_decay << endl;
    fidFile << "K_spont_ass=" << (double)K_spont_ass << endl;
    fidFile << "K_spont_diss=" << (double)K_spont_diss << endl;
    fidFile << "moleculeDecay_KineticConstant=" << (double)moleculeDecay_KineticConstant << endl << endl;

    fidFile << "# (0 or 0.5) if set to 0.5 the speed of molecules goes with the inverse of the square of the length" << endl;
    fidFile << "diffusion_contribute=" << (double)diffusion_contribute << endl << endl;

    fidFile << "# Solubility Threshold (The threshold set the centre of the logistic curve, 0 to exclude precipitation)" << endl;
    fidFile << "solubility_threshold=" << solubility_threshold << endl << endl;

    fidFile << "# Overall influx (moles per second) and maximum length of the species passing the filter in the outflow process" << endl;
    fidFile << "influx_rate=" << (double)influx_rate << endl;
    fidFile << "maxLOut=" << maxLOut << endl << endl;


    fidFile.close();
    cout << "done." << endl;

    if(debugLevel == FINDERRORDURINGRUNTIME) cout << "environment::saveConfigurationFileSTD end" << endl;

    return true;
}

/**
 Save influx structure on file, standard C++
 @param bool saveInfluxStructure(QString tmpStoringPath);
 @version 1.0
 @date 2013-047-03
*/
bool environment::saveInfluxStructureSTD(string tmpStoringPath)
{
    if(debugLevel == FINDERRORDURINGRUNTIME) cout << "environment::saveInfluxStructureSTD start" << endl;
    cout << "Saving Influx file...";

    string strInfluxFile = tmpStoringPath + "/_acsinflux.csv";

    ofstream fidFile;

    try{
        fidFile.open(strInfluxFile.c_str(), ios::out | ios::app);
    }
    catch(exception&e)
    {
        cerr << "exceptioncaught:" << e.what() << endl;
        ExitWithError("error in method saveInfluxStructureSTD","exceptionerrorthrown");
    }

    acs_int influxSpeciesID = 0;
    for(vector<acs_int>::iterator tmpAllInfluxIter = nutrientsForInflux.begin(); tmpAllInfluxIter != nutrientsForInflux.end(); tmpAllInfluxIter++)
    {
        fidFile << *tmpAllInfluxIter << "\t" << (double)nutrientsProb2BeSelected.at(influxSpeciesID) << endl;
        influxSpeciesID++;
    }
    fidFile.close();
    cout << "done." << endl;

    if(debugLevel == FINDERRORDURINGRUNTIME) cout << "environment::saveInfluxStructureSTD end" << endl;
    return true;

}

/**
 Save Energetic Boolean Function on file, standard C++
 @param bool saveInfluxStructure(QString tmpStoringPath);
 @version 1.0
 @date 2013-047-03
*/
bool environment::saveNrgBoolFncStructureSTD(string tmpStoringPath)
{
    if(debugLevel == FINDERRORDURINGRUNTIME) cout << "environment::saveNrgBoolFncStructureSTD start" << endl;
    cout << "Saving Influx file...";

    string strNrgFile = tmpStoringPath + "/_acsnrgbooleanfunctions.csv";

    ofstream fidFile;

    try{
        fidFile.open(strNrgFile.c_str(), ios::out | ios::app);
    }
    catch(exception&e)
    {
        cerr << "exceptioncaught:" << e.what() << endl;
        ExitWithError("error in method saveNrgBoolFncStructureSTD","exceptionerrorthrown");
    }

    acs_int influxNrgBoolID = 0;
    for(vector<acs_int>::iterator tmpAllNrgBoolIter = nrgBooleanFunctions.begin(); tmpAllNrgBoolIter != nrgBooleanFunctions.end(); tmpAllNrgBoolIter++)
    {
        fidFile << *tmpAllNrgBoolIter << "\t" << (double)nrgBoolFncsProb2BeSelected.at(influxNrgBoolID) << endl;
        influxNrgBoolID++;
    }
    fidFile.close();
    cout << "done." << endl;

    if(debugLevel == FINDERRORDURINGRUNTIME) cout << "environment::saveNrgBoolFncStructureSTD end" << endl;
    return true;

}

/**
 Save the species structures in a file named species_[currentSims]_[currentStep].csv. This is file is equal to the acs_species.csv input file. C++ standard
 The file is saved in the directory indicated as a second parameter in the run command
 @version 1.0
 @date 2013/07/03
 */
bool environment::saveSpeciesStructureSTD(acs_int tmpCurrentGen, acs_int tmpCurrentSim, acs_int tmpCurrentStep, string tmpStoringPath)
{
    if(debugLevel == FINDERRORDURINGRUNTIME) cout << "environment::saveSpeciesStructure start" << endl;

    if(debugLevel >= SMALL_DEBUG)
            cout << "Saving species structure...";
    stringstream strCurrenGen;
    stringstream strCurrentSim;
    stringstream strCurrentStep;
    string strZeroGenBefore = zeroBeforeStringNumberSTD (nGEN, tmpCurrentGen);
    string strZeroSimBefore = zeroBeforeStringNumberSTD(pow(double(nSIM), double(tmpCurrentGen)), tmpCurrentSim);
    string strZeroStepBefore = zeroBeforeStringNumberSTD(nReactions, tmpCurrentStep);

    strCurrenGen << tmpCurrentGen;
    strCurrentSim << tmpCurrentSim;
    strCurrentStep << tmpCurrentStep;

    string strFileSpeciesStructure = tmpStoringPath + "/species_" + strZeroGenBefore + strCurrenGen.str() + "_" +
                                    strZeroSimBefore + strCurrentSim.str() + "_" +
                                    strZeroStepBefore + strCurrentStep.str() + ".csv";

    ofstream fidFile;
    try{
        fidFile.open(strFileSpeciesStructure.c_str(), ios::out | ios::app);
    }
    catch(exception&e)
    {
        cerr << "exceptioncaught:" << e.what() << endl;
        ExitWithError("error in method saveSpeciesStructureSTD","exceptionerrorthrown");
    }

    for(vector<species>::iterator tmpAllSpeciesIter = allSpecies.begin(); tmpAllSpeciesIter != allSpecies.end(); tmpAllSpeciesIter++)
    {
        fidFile << (acs_longInt)tmpAllSpeciesIter->getID() << "\t"
        << tmpAllSpeciesIter->getSequence().c_str() << "\t"
        //<< (acs_longInt)tmpAllSpeciesIter->getAmount() << "\t"
        << (double)tmpAllSpeciesIter->getConcentration() << "\t"
        << (double)tmpAllSpeciesIter->getDiffusionEnh() << "\t"
        << (acs_int)tmpAllSpeciesIter->getSolubility() << "\t"
        << (double)tmpAllSpeciesIter->getComplexDegEnh() << "\t"
        << tmpAllSpeciesIter->getComplexCutPnt() << "\t"
        << tmpAllSpeciesIter->getEvaluated() << "\t"
        << (double)tmpAllSpeciesIter->getAge() << "\t"
        << tmpAllSpeciesIter->getReborns() << "\t"
        << tmpAllSpeciesIter->getCatalyst_ID()  << "\t"
        << tmpAllSpeciesIter->getSubstrate_ID() <<	"\t"
        << (double)tmpAllSpeciesIter->getK_phospho() << "\t"
        << (double)tmpAllSpeciesIter->getLoadedConcentration(volume) << "\t"
        << tmpAllSpeciesIter->getConcentrationFixed() << endl;
    }
    fidFile.close();

    if(debugLevel >= SMALL_DEBUG)
        cout << "done." << endl;

    if(debugLevel == FINDERRORDURINGRUNTIME) cout << "environment::saveSpeciesStructure end" << endl;

    return true;
}

/**
 Save the reactions structures in a file named reactions_[currentSims]_[currentStep].csv. This is file is equal to the acs_reactions.csv input file.
 The file is saved in the directory indicated as a second parameter in the run command
 @version 1.0
 */
bool environment::saveReactionsStructureSTD(acs_int tmpCurrentGen, acs_int tmpCurrentSim, acs_int tmpCurrentStep, string tmpStoringPath)
{
    if(debugLevel == FINDERRORDURINGRUNTIME) cout << "environment::saveReactionsStructureSTD start" << endl;

    if(debugLevel >= SMALL_DEBUG)
        cout << "Saving reactions structure...";
    stringstream  strCurrenGen;
    stringstream strCurrentSim;
    stringstream strCurrentStep;
    string strZeroGenBefore = zeroBeforeStringNumberSTD(nGEN, tmpCurrentGen);
    string strZeroSimBefore = zeroBeforeStringNumberSTD(pow(double(nSIM), double(tmpCurrentGen)), tmpCurrentSim);
    string strZeroStepBefore = zeroBeforeStringNumberSTD(nReactions, tmpCurrentStep);
    strCurrenGen << tmpCurrentGen;
    strCurrentSim << tmpCurrentSim;
    strCurrentStep << tmpCurrentStep;
    string strFileReactionsStructure = tmpStoringPath + "/reactions_" + strZeroGenBefore + strCurrenGen.str() + "_" +
            strZeroSimBefore + strCurrentSim.str() + "_" +
            strZeroStepBefore + strCurrentStep.str() + ".csv";

    ofstream fidFile;
    try{
        fidFile.open(strFileReactionsStructure.c_str(), ios::out | ios::app);
    }
    catch(exception&e)
    {
        cerr << "exceptioncaught:" << e.what() << endl;
        ExitWithError("error in method saveReactionsStructureSTD","exceptionerrorthrown");
    }

    for(vector<reactions>::iterator tmpAllReactionsIter = allReactions.begin(); tmpAllReactionsIter != allReactions.end(); tmpAllReactionsIter++)
    {
        fidFile << tmpAllReactionsIter->getID() << "\t"
        << tmpAllReactionsIter->getType() << "\t"
        << tmpAllReactionsIter->getSpecies_I() << "\t"
        << tmpAllReactionsIter->getSpecies_II() << "\t"
        << tmpAllReactionsIter->getSpecies_III() << "\t"
        << tmpAllReactionsIter->getEvents() << "\t"
        << tmpAllReactionsIter->getEnergyType() << "\t"
        << (acs_double)tmpAllReactionsIter->getKspont() << endl;
    }
    fidFile.close();
    if(debugLevel >= SMALL_DEBUG)
        cout << "done." << endl;

    if(debugLevel == FINDERRORDURINGRUNTIME) cout << "environment::saveReactionsStructureSTD end" << endl;

    return true;
}

/**
 Save the catalysis structures in a file named catalysis_[currentSims]_[currentStep].csv. This is file is equal to the acs_catalysis.csv input file. Standard C++
 The file is saved in the directory indicated as a second parameter in the run command
 @version 1.0
 */
bool environment::saveCatalysisStructureSTD(acs_int tmpCurrentGen, acs_int tmpCurrentSim, acs_int tmpCurrentStep, string tmpStoringPath)
{
    if(debugLevel == FINDERRORDURINGRUNTIME) cout << "environment::saveCatalysisStructureSTD start" << endl;

    if(debugLevel >= SMALL_DEBUG)
        cout << "Saving catalysis structure...";
    stringstream strCurrenGen;
    stringstream strCurrentSim;
    stringstream strCurrentStep;
    string strZeroGenBefore = zeroBeforeStringNumberSTD(nGEN, tmpCurrentGen);
    string strZeroSimBefore = zeroBeforeStringNumberSTD(pow(double(nSIM), double(tmpCurrentGen)), tmpCurrentSim);
    string strZeroStepBefore = zeroBeforeStringNumberSTD(nReactions, tmpCurrentStep);

    strCurrenGen << tmpCurrentGen;
    strCurrentSim << tmpCurrentSim;
    strCurrentStep << tmpCurrentStep;

    string strFileCatalysisStructure = tmpStoringPath + "/catalysis_" + strZeroGenBefore + strCurrenGen.str() + "_" +
            strZeroSimBefore + strCurrentSim.str() + "_" +
            strZeroStepBefore + strCurrentStep.str() + ".csv";

    ofstream fidFile;
    try{
        fidFile.open(strFileCatalysisStructure.c_str(), ios::out | ios::app);
    }
    catch(exception&e)
    {
        cerr << "exceptioncaught:" << e.what() << endl;
        ExitWithError("error in method saveCatalysisStructureSTD","exceptionerrorthrown");
    }

    for(vector<catalysis>::iterator tmpAllCatalysisIter = allCatalysis.begin(); tmpAllCatalysisIter != allCatalysis.end(); tmpAllCatalysisIter++)
    {
        fidFile << tmpAllCatalysisIter->getCatId() << "\t"
        << tmpAllCatalysisIter->getCat() << "\t"
        << tmpAllCatalysisIter->getReactionID() << "\t"
        << tmpAllCatalysisIter->getTotAmount() << "\t"
        << (double)tmpAllCatalysisIter->getKass() << "\t"
        << (double)tmpAllCatalysisIter->getKdiss() << "\t"
        << (double)tmpAllCatalysisIter->getK_cpx() << "\t"
        << tmpAllCatalysisIter->getCpxTarget() << endl;
    }
    fidFile.close();
    if(debugLevel >= SMALL_DEBUG)
        cout << "done." << endl;

    if(debugLevel == FINDERRORDURINGRUNTIME) cout << "environment::saveCatalysisStructureSTD end" << endl;

    return true;
}

/**
 Save the reactions times in a file named times_[currentSim].csv. Standard C++
 The file is saved in the directory indicated as a second parameter in the run command
 @version 1.0
 @date 2013/07/03
 */
bool environment::saveTimesSTD(acs_int tmpCurrentGen, acs_int tmpCurrentSim, acs_int tmpCurrentStep, string tmpStoringPath)
{
    if(debugLevel == FINDERRORDURINGRUNTIME) cout << "environment::saveTimesSTD start" << endl;
    if(debugLevel >= SMALL_DEBUG)
        cout << "\t|- Saving Times to file...";

    try{

        stringstream strCurrentGen;
        stringstream strCurrentSim;
        string strZeroGenBefore = zeroBeforeStringNumberSTD(nGEN, tmpCurrentGen);
        string strZeroSimBefore = zeroBeforeStringNumberSTD(pow(double(nSIM), double(tmpCurrentGen)), tmpCurrentSim);

        strCurrentGen << tmpCurrentGen;
        strCurrentSim << tmpCurrentSim;

        string strFileTimesStructure = tmpStoringPath + "/times_" + strZeroGenBefore + strCurrentGen.str() + "_" +
                strZeroSimBefore + strCurrentSim.str() + ".csv";

        ofstream fidFile;
        try{
            fidFile.open(strFileTimesStructure.c_str(), ios::out | ios::app);
        }
        catch(exception&e)
        {
            cerr << "exceptioncaught:" << e.what() << endl;
            ExitWithError("error in method saveTimesSTD","exceptionerrorthrown");
        }

        if(COPYOFallGillespieScores.size() > 0)
        {
            fidFile << tmpCurrentStep << "\t"
                << (double)reactionsTime.at(tmpCurrentStep-1) << "\t"
                << gillespieReactionsSelected.at(tmpCurrentStep-1) << "\t"
                << COPYOFallGillespieScores.at(gillespieReactionsSelected.at(tmpCurrentStep-1)).getIdReactionType() << "\t"
                << COPYOFallGillespieScores.size() << "\t"
                << (double)allTimes.at(tmpCurrentStep-1) << "\t"
                << getTotalNumberOfSpecies() << "\t"
                << getTotalNumberOfMolecules() << "\t"
                << getTotalNumberOfComplexSpecies() << "\t"
                << getTotalNumberOfComplexes() << "\t"
                << getTotalNumberOfMonomers() << "\t"
                << (double)gillespiePartialTimes.at(tmpCurrentStep-1) << "\t"
                << (double)performReactionPartialTimes.at(tmpCurrentStep-1) << "\t"
                << (double)remainingProcessesPartialTimes.at(tmpCurrentStep-1) << "\t"
                << (double)ratioBetweenNewGillTotGill << endl;
        }else{
            fidFile << tmpCurrentStep << "\t"
                << (double)reactionsTime.at(tmpCurrentStep-1) << "\t"
                << gillespieReactionsSelected.at(tmpCurrentStep-1) << "\t"
                << 0 << "\t"
                << 0 << "\t"
                << 0 << "\t"
                << getTotalNumberOfSpecies() << "\t"
                << getTotalNumberOfMolecules() << "\t"
                << getTotalNumberOfComplexSpecies() << "\t"
                << getTotalNumberOfComplexes() << "\t"
                << getTotalNumberOfMonomers() << "\t"
                << 0 << "\t"
                << 0 << "\t"
                << 0 << "\t"
                << 0 << endl;
        }
        fidFile.close();
        if(debugLevel >= SMALL_DEBUG)
            cout << "OK" << endl;

        COPYOFallGillespieScores.clear();

    }
    catch(exception&e)
    {
        cerr << "exceptioncaught:" << e.what() << endl;
        ExitWithError("error in method saveTimes","exceptionerrorthrown");
    }

    if(debugLevel == FINDERRORDURINGRUNTIME) cout << "environment::saveTimesSTD end" << endl;

    return true;
}

/**
 Save the reactions parameters in a file named reactions_parameters_[currentSim].csv.
 The file is saved in the directory indicated as a second parameter in the run command
 @version 1.0
 @date 2013/07/03
 */
bool environment::saveReactionsParametersSTD(acs_int tmp__CurrentGen, acs_int tmp__CurrentSim, acs_int tmp__CurrentStep, string tmp__StoringPath,
                                          acs_int tmpRctType, acs_longInt tmpCat, acs_longInt tmpMol_I, acs_longInt tmpMol_II, acs_longInt tmpMol_III)
{
    if(debugLevel == FINDERRORDURINGRUNTIME) cout << "environment::saveReactionsParametersSTD start" << endl;
    if(debugLevel >= HIGH_DEBUG)
            cout << "\t|- Saving reaction parameters to file...";

    stringstream strCurrentGen;
    stringstream strCurrentSim;
    string strZeroGenBefore = zeroBeforeStringNumberSTD(nGEN, tmp__CurrentGen);
    string strZeroSimBefore = zeroBeforeStringNumberSTD(pow(double(nSIM), double(tmp__CurrentGen)), tmp__CurrentSim);

    strCurrentGen << tmp__CurrentGen;
    strCurrentSim << tmp__CurrentSim;

    string strFileReactionsParameters = tmp__StoringPath + "/reactions_parameters_" + strZeroGenBefore + strCurrentGen.str() + "_" +
            strZeroSimBefore + strCurrentSim.str() + ".csv";

    ofstream fidFile;
    try{
        fidFile.open(strFileReactionsParameters.c_str(), ios::out | ios::app);
    }
    catch(exception&e)
    {
        cerr << "exceptioncaught:" << e.what() << endl;
        ExitWithError("error in method saveTimesSTD","exceptionerrorthrown");
    }

    fidFile << tmp__CurrentStep << "\t"
    << (double)reactionsTime.at(tmp__CurrentStep-1) << "\t"
    << tmpRctType << "\t"
    << tmpCat << "\t"
    << tmpMol_I << "\t"
    << tmpMol_II << "\t"
    << tmpMol_III << "\t"
    << overallLoadedMolsCounter << "\t"
    << (double)(overallLoadedMolsCounter / (AVO * volume)) << "\t"
    << (double)gillespieMean << "\t"
    << (double)gillespieSD << "\t"
    << (double)gillespieEntropy << "\t"
    << (double)ratioBetweenNewGillTotGill << "\t"
    << (double)ratioBetweenReverseAndTotalScore << endl;
    fidFile.close();
    if(debugLevel >= HIGH_DEBUG)
        cout << "OK" << endl;

    COPYOFallGillespieScores.clear();

    if(debugLevel == FINDERRORDURINGRUNTIME) cout << "environment::saveReactionsParametersSTD end" << endl;

    return true;
}

/**
 Save living species in a file named living_species_[currentSim].csv. Standard C++
 The file is saved in the directory indicated as a second parameter in the run command
 @version 1.0
 @date 2013/07/03
 */
bool environment::saveLivingSpeciesIDSTD(acs_int tmp__CurrentGen, acs_int tmp__CurrentSim, acs_int tmp__CurrentStep, string tmp__StoringPath)
{
    if(debugLevel == FINDERRORDURINGRUNTIME) cout << "environment::saveLivingSpeciesIDSTD start" << endl;
    try{

        if(debugLevel == FINDERRORDURINGRUNTIME) cout << "environment::saveLivingSpeciesIDSTD start" << endl;
        if(debugLevel >= SMALL_DEBUG)
            cout << "\t|- Saving saveLivingSpeciesIDSTD to file...";

        stringstream strCurrentGen;
        stringstream strCurrentSim;
        string strZeroGenBefore = zeroBeforeStringNumberSTD(nGEN, tmp__CurrentGen);
        string strZeroSimBefore = zeroBeforeStringNumberSTD(pow(double(nSIM), double(tmp__CurrentGen)), tmp__CurrentSim);

        strCurrentGen << tmp__CurrentGen;
        strCurrentSim << tmp__CurrentSim;

        string strFilelivingSpecies = tmp__StoringPath + "/livingSpecies_" + strZeroGenBefore + strCurrentGen.str() + "_" +
                strZeroSimBefore + strCurrentSim.str() + ".csv";

        ofstream fidFile;
        try{
            fidFile.open(strFilelivingSpecies.c_str(), ios::out | ios::app);
        }
        catch(exception&e)
        {
            cerr << "exceptioncaught:" << e.what() << endl;
            ExitWithError("error in method saveTimesSTD","exceptionerrorthrown");
        }

        fidFile << tmp__CurrentStep << "\t" << (double)reactionsTime.at(tmp__CurrentStep-1);
        for(acs_longInt i = 0; i < (acs_longInt)allSpecies.size(); i++)
        {
            if((allSpecies.at(i).getAmount() > 0) & (allSpecies.at(i).getComplexCutPnt() == 0))
                fidFile << "\t" << allSpecies.at(i).getID();
        }
        fidFile << endl;
        fidFile.close();
        if(debugLevel >= SMALL_DEBUG)
            cout << "OK" << endl;

        COPYOFallGillespieScores.clear();

        if(debugLevel == FINDERRORDURINGRUNTIME) cout << "environment::saveLivingSpeciesIDSTD end" << endl;

    }
    catch(exception&e)
    {
        cerr << "exceptioncaught:" << e.what() << endl;
        ExitWithError("error in saveLivingSpeciesIDSTD method","exceptionerrorthrown");
    }

    return true;
}

/**
 Save living species total AMOUNT in a file named livingAmount_[CurrentGen]_[currentSim].csv.
 The file is saved in the directory indicated as a second parameter in the run command - Standard C++
 @version 1.0
 @date 2013/07/03
 */
bool environment::saveLivingSpeciesAmountSTD(acs_int tmp__CurrentGen, acs_int tmp__CurrentSim, string tmp__StoringPath)
{
    try{
        if(debugLevel == FINDERRORDURINGRUNTIME) cout << "environment::saveLivingSpeciesAmountSTD start" << endl;
        if(debugLevel >= SMALL_DEBUG)
            cout << "\t|- Saving reaction parameters to file...";

        stringstream strCurrentGen;
        stringstream strCurrentSim;
        string strZeroGenBefore = zeroBeforeStringNumberSTD(nGEN, tmp__CurrentGen);
        string strZeroSimBefore = zeroBeforeStringNumberSTD(pow(double(nSIM), double(tmp__CurrentGen)), tmp__CurrentSim);
        strCurrentGen << tmp__CurrentGen;
        strCurrentSim << tmp__CurrentSim;

        string strFilelivingSpecies = tmp__StoringPath + "/livingAmount_" + strZeroGenBefore + strCurrentGen.str() + "_" +
                strZeroSimBefore + strCurrentSim.str() + ".csv";

        ofstream fidFile;
        try{
            fidFile.open(strFilelivingSpecies.c_str(), ios::out | ios::app);
        }
        catch(exception&e)
        {
            cerr << "exceptioncaught:" << e.what() << endl;
            ExitWithError("error in method saveTimesSTD","exceptionerrorthrown");
        }

        acs_int firstSpeciesControl = 0;
        for(acs_longInt i = 0; i < (acs_longInt)allSpecies.size(); i++)
        {
            if((allSpecies.at(i).getAmount() > 0) & (allSpecies.at(i).getComplexCutPnt() == 0))
            {
                if(firstSpeciesControl != 0)
                {
                    fidFile << "\t" << allSpecies.at(i).getAmount();
                }else{
                    fidFile << allSpecies.at(i).getAmount();
                    firstSpeciesControl++;
                }
            }
        }
        fidFile << endl;
        fidFile.close();
        if(debugLevel >= SMALL_DEBUG)
            cout << "OK" << endl;

        COPYOFallGillespieScores.clear();

        if(debugLevel == FINDERRORDURINGRUNTIME) cout << "environment::saveLivingSpeciesAmountSTD end" << endl;

    }
    catch(exception&e)
    {
        cerr << "exceptioncaught:" << e.what() << endl;
        ExitWithError("error in saveLivingSpeciesAmount method","exceptionerrorthrown");
    }

    return true;
}

/**
 Save living species total CONCENTRATION in a file named livingAmount_[CurrentGen]_[currentSim].csv.
 The file is saved in the directory indicated as a second parameter in the run command
 @version 1.0
 @date 2013/07/03
 */
bool environment::saveLivingSpeciesConcentrationSTD(acs_int tmp__CurrentGen, acs_int tmp__CurrentSim, string tmp__StoringPath)
{
    try{
        if(debugLevel == FINDERRORDURINGRUNTIME) cout << "environment::saveLivingSpeciesConcentrationSTD start" << endl;
        if(debugLevel >= SMALL_DEBUG)
            cout << "\t|- Saving reaction parameters to file...";

        stringstream strCurrentGen;
        stringstream strCurrentSim;
        string strZeroGenBefore = zeroBeforeStringNumberSTD(nGEN, tmp__CurrentGen);
        string strZeroSimBefore = zeroBeforeStringNumberSTD(pow(double(nSIM), double(tmp__CurrentGen)), tmp__CurrentSim);
        strCurrentGen << tmp__CurrentGen;
        strCurrentSim << tmp__CurrentSim;
        string strFilelivingSpecies = tmp__StoringPath + "/livingConcentration_" + strZeroGenBefore + strCurrentGen.str() + "_" +
                strZeroSimBefore + strCurrentSim.str() + ".csv";

        ofstream fidFile;
        try{
            fidFile.open(strFilelivingSpecies.c_str(), ios::out | ios::app);
        }
        catch(exception&e)
        {
            cerr << "exceptioncaught:" << e.what() << endl;
            ExitWithError("error in method saveTimesSTD","exceptionerrorthrown");
        }
        acs_int firstSpeciesControl = 0;
        for(acs_longInt i = 0; i < (acs_longInt)allSpecies.size(); i++)
        {
            if((allSpecies.at(i).getAmount() > 0) & (allSpecies.at(i).getComplexCutPnt() == 0))
            {
                if(firstSpeciesControl != 0)
                {
                    fidFile << "\t" << (double)allSpecies.at(i).getConcentration();
                }else{
                    fidFile << (double)allSpecies.at(i).getConcentration();
                    firstSpeciesControl++;
                }
            }
        }
        fidFile << endl;
        fidFile.close();
        if(debugLevel >= SMALL_DEBUG)
            cout << "OK" << endl;

        COPYOFallGillespieScores.clear();

        if(debugLevel == FINDERRORDURINGRUNTIME) cout << "environment::saveLivingSpeciesConcentrationSTD end" << endl;

    }
    catch(exception&e)
    {
        cerr << "exceptioncaught:" << e.what() << endl;
        ExitWithError("error in saveLivingSpeciesConcentration method","exceptionerrorthrown");
    }

    return true;
}

///////////////////////////////////
// STATISTIC FUNTIONS
///////////////////////////////////

bool environment::devStd()
{
    if(debugLevel == FINDERRORDURINGRUNTIME) cout << "environment::devStd start" << endl;
    bool devStdFlag = true;
    acs_double tempDevStd = 0;

    for(vector<gillespie>::iterator tmpScoreValue = allGillespieScores.begin(); tmpScoreValue != allGillespieScores.end(); tmpScoreValue++)
    {
        tempDevStd += pow((double)tmpScoreValue->getScore() - gillespieMean,2);
    }

    gillespieSD = sqrt((double)tempDevStd / (double)allGillespieScores.size());

    if(debugLevel == FINDERRORDURINGRUNTIME) cout << "environment::devStd end" << endl;

    return devStdFlag;
}

bool environment::entropy()
{
    if(debugLevel == FINDERRORDURINGRUNTIME) cout << "environment::entropy start" << endl;
    bool entropyFlag = true;
    acs_double tempEntropy = 0;
    acs_double normValue;

    for(vector<gillespie>::iterator tmpScoreValue = allGillespieScores.begin(); tmpScoreValue != allGillespieScores.end(); tmpScoreValue++)
    {
        normValue = tmpScoreValue->getScore() / gillespieCumulativeStepScoreList.back();
        tempEntropy += -normValue * log(normValue);
    }

    gillespieEntropy = tempEntropy;

    if(debugLevel == FINDERRORDURINGRUNTIME) cout << "environment::entropy end" << endl;

    return entropyFlag;
}

