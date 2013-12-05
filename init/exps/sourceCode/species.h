/** class species
 *	\brief This class contains declarations of the species class
 *
 *  @authors Alessandro Filisetti
 *  @version 1.1 questa modifica è di prova per subversion
 *
 *  Created by Alessandro Filisetti on 19/02/09.
 *  Copyright 2009 European Centre for Living Technology. All rights reserved.
 *  Test paxelito SVN
 */

#ifndef SPECIES
#define SPECIES

#include "acs_headers.h"
#include "commonFunctions.h"

class species
{
private:
    acs_longInt id;                     /**< Species identificator */
    string sequence;                    /**< Species sequence */
	acs_longInt amount;					/**< Total Number of molecules (charged + unchurged) belonging to this species */
	acs_longInt chargedMols;			/**< Number of charged molecules*/
	acs_double concentration;			/**< Concentration of the species */
	acs_double age;						/**< Species age (computed from the last creation) */
	acs_int reborns;					/**< Number of times that the species became greater than zero */
	acs_double diffusionEnh;			/**< Diffusion initialization parameter */
	acs_int soluble;					/**< Solubility parameter Threshold: 1 soluble - 0 precipitated (this species will not be avalaible for reactions) */
	acs_double complexDegradationEnh;	/**< Complex degradation rate (If the species is not a complex this value is set to 0) */
	acs_int complexCuttingPoint;		/**< Binding point in case of complex (3 means for example that ABBAA could be formed by a catalyst AAB and a substrate AA) */
	acs_int evaluated;					/**< Boolean variable that indicates the state of the species (1: already evaluated, 0: virtual) */
	acs_longInt catalyst_ID;			/**< If the species is a complex this is the ID of the catalyst */
	acs_longInt substrate_ID;			/**< If the species is a complex this is the ID of the substrate */
	acs_double K_phospho;				/**< Energy Phosphorilation kinetic constant associated with each species */
	acs_int energizable;				/**< This parameter indicates whether the species is energizable or not*/
	bool concentrationFixed;			/**< This parameter indicates whether or not the species concentration may change*/
    acs_double firstConcentration;      /**< Concentration at the beginning of the simulation, > 0 for the initial species, 0 otherwise*/
    acs_int lastSpeciesEvaluated;		/**< ID of the last species evaluated by this species for possible reactions*/
    vector<acs_longInt> secondSubstrates;/**< If the species in a complex this is a list of the second substrates with which perform the final condensation*/
    vector<acs_longInt> catalysisIfCpx; /**< Catalysis ID list if complex*/
    vector<acs_double> kCond;			/**< List parallel to secondSubstrates, contain a copy of the k condensation*/
    vector<acs_longInt> gillespieEngagement; /**< List containing the Gillespie's ID in which the species is involved */
	
public:
	//!< New species constructor (IN AMOUNT)
	species(); //Default species constructor, it is used to create the energy carrier
	species(acs_longInt tmpID, string tmpSequence, acs_longInt tmpAmount, acs_double tmpDiffusionEnh, 
			acs_int tmpSoluble, acs_double tmpComplexDegEnh, acs_int tmpComplexCuttingPoint,
			acs_int tmpEvalueted, acs_double tmpVolume, acs_double tmpK_phospho, acs_int tmpEnergizable, 
			acs_double tmpInflux_rate, acs_int tmpMaxLOut);
	
	//!< New species constructor (IN CONCENTRATION)
	species(acs_longInt tmpID, string tmpSequence, acs_double tmpConcentration, acs_double tmpDiffusionEnh, 
			acs_int tmpSoluble, acs_double tmpComplexDegEnh, acs_int tmpComplexCuttingPoint,
			acs_int tmpEvalueted, acs_double tmpVolume, acs_double tmpK_phospho, acs_int tmpEnergizable, 
			acs_double tmpInflux_rate, acs_int tmpMaxLOut);
	
	//!< New species constructor in case of species structure file upload (IN AMOUNT)
	species(acs_longInt tmpID, string tmpSequence, acs_longInt tmpAmount, acs_double tmpDiffusionEnh, 
			acs_int tmpSoluble, acs_double tmpComplexDegEnh, acs_int tmpComplexCuttingPoint,
			acs_int tmpEvalueted, acs_double tmpAge, acs_int tmpReborns, acs_double tmpVolume, 
			acs_longInt tmpNotUsedCatID, acs_longInt tmpNotUsedSubID, acs_double tmpK_phospho, acs_int tmpEnergizable, 
			acs_double tmpInflux_rate, acs_int tmpMaxLOut);
	
	//!< New species constructor in case of species structure file upload (IN CONCENTRATION)
	species(acs_longInt tmpID, string tmpSequence, acs_double tmpConcentration, acs_double tmpDiffusionEnh, 
			acs_int tmpSoluble, acs_double tmpComplexDegEnh, acs_int tmpComplexCuttingPoint,
			acs_int tmpEvalueted, acs_double tmpAge, acs_int tmpReborns, acs_double tmpVolume, 
			acs_longInt tmpNotUsedCatID, acs_longInt tmpNotUsedSubID, acs_double tmpK_phospho, acs_double tmpKLoadConc, 
			acs_int tmpEnergizable, acs_double tmpInflux_rate, acs_int tmpMaxLOut);	
	
	//!< New random species constructor
	species(acs_longInt tmpID, string tmpSequence, acs_longInt tmpAmount, acs_double tmpDiffusionEnh, 
			acs_int tmpSoluble, acs_double tmpComplexProb, acs_double tmpMaxComplexDegKinetic,
			MTRand& tmp_RandomGenerator, acs_double tmpVolume, acs_double tmpK_phospho, acs_int tmpEnergizable);
	
	//!< new Complex species constructor
	species(acs_longInt tmpID, string tmpSequence, acs_double tmpDiffusionEnh, 
			acs_int tmpSoluble, acs_double tmpMaxComplexDegKinetic, 
			acs_int tmpCuttingPoint, MTRand& tmp_RandomGenerator, acs_longInt tmpCatalyst_ID, 
			acs_longInt tmpSubstrate_ID, acs_double tmpVolume, acs_double tmpK_phospho, acs_int tmpEnergizable);
	~species(){}
	
	// GETTING FUNCTIONS
	acs_longInt getID()const{return id;}
	string getSequence()const{ return sequence;}
	acs_int getSequenceLength()const{return (unsigned)sequence.length();}
	acs_longInt getAmount()const{return amount;}
	acs_longInt getNOTchargeMols()const{return amount - chargedMols;}
	acs_longInt getChargeMols()const{return chargedMols;}
	acs_double getConcentration()const{return concentration;}
	acs_double getLoadedConcentration(acs_double tmpVolume);
	acs_double getAge()const{return age;}
	acs_int getReborns()const{return reborns;} 
	acs_double getDiffusionEnh()const{return diffusionEnh;}
	acs_int getSolubility()const{return soluble;}
	acs_double getComplexDegEnh()const{return complexDegradationEnh;}
	acs_int getComplexCutPnt()const{return complexCuttingPoint;}
	acs_int getEvaluated()const{return evaluated;}
	acs_longInt getCatalyst_ID()const{return catalyst_ID;}
	acs_longInt getSubstrate_ID()const{return substrate_ID;}
	acs_double getK_phospho()const{return K_phospho;}
	acs_int getEnergizable()const{return energizable;}
	bool getConcentrationFixed()const{return concentrationFixed;}
    acs_double getFirstConcentration()const{return firstConcentration;}
    acs_int getLastSpeciesEvaluated()const{return lastSpeciesEvaluated;}
    vector<acs_longInt> getSecSubList()const{return secondSubstrates;}
    vector<acs_double> getSec_k_SubList()const{return kCond;}
    vector<acs_longInt> getCatalysisIfCpx()const{return catalysisIfCpx;}
    acs_int getSecSubListSize()const{return secondSubstrates.size();}
    acs_longInt getSecSubListID(acs_int tmpID)const{return secondSubstrates.at(tmpID);}
    acs_double getSec_k_SubListID(acs_int tmpID)const{return kCond.at(tmpID);}
    acs_longInt getCatalysisIfCpxID(acs_int tmpID)const{return catalysisIfCpx.at(tmpID);}
	
	// FUNCTIONAL FUNCTIONS
	// Total quantity update functions (concentration is recomputed)
    void increment(acs_double tmpVolume);
	void specificIncrement(acs_int tmpIncrement, acs_double tmpVolume){amount = amount + tmpIncrement; numToConc(tmpVolume);}
	void setAmount(acs_int tmpAmount, acs_double tmpVolume){amount = tmpAmount; numToConc(tmpVolume);}
	void setConcentration(acs_double tmpConc, acs_double tmpVolume){concentration = tmpConc; concToNum(tmpVolume);}
    //TR void decrement(acs_double tmpVolume){amount--; numToConc(tmpVolume);}
    void decrement(acs_double tmpVolume);
	bool setChargeMols(acs_int tmpMolsToCharge);
	bool setSpecificChargeMols(acs_int tmpMolsToCharge);
	bool chargeMol();
	bool unchargeMol();
	// set new values
	void setEvaluated(){evaluated = 1;}
	void setDiffusion(acs_double tmpDiff){diffusionEnh = tmpDiff;}
	void setSolubility(acs_int tmpSol){soluble = tmpSol;}
	void setKphospho(acs_double tmpKphospho){K_phospho = tmpKphospho;}
	void setNewAge(acs_double tmpLastTimeInterval){age += tmpLastTimeInterval;}
	void rebornsIncrement(){reborns++;}
	// Concentration and amount update
	void concToNum(acs_double tmpVolume){amount = acsround(AVO * concentration * tmpVolume);}
	void numToConc(acs_double tmpVolume){concentration = amount / (AVO * tmpVolume);}
    // Reset Functions
    void resetAge(){age=0;}
    void resetReborns(){reborns=0;}
    void resetToInitConc(acs_double tmpVolume){concentration=firstConcentration; concToNum(tmpVolume);}
    void setLastSpeciesEvaluated(acs_int tmpID){lastSpeciesEvaluated = tmpID;}
    void insertSecSub(acs_longInt tmpID, acs_double tmpK, acs_longInt tmpCat);
    void insertGillID(acs_longInt tmpID){gillespieEngagement.push_back(tmpID);}
    acs_longInt getGillIDpos(acs_longInt tmpID)const{return gillespieEngagement.at(tmpID);}
    void showGillEngagement();
};

#endif
