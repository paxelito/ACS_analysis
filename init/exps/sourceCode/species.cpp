/** class species
 *	\brief This class containing the declaration of the species
 *
 *  @authors Alessandro Filisetti
 *  @version 0.1
 *
 *  Created by Alessandro Filisetti on 19/02/09.
 *  Copyright 2009 European Centre for Living Technology. All rights reserved.
 *
 */

#include "species.h"

/** \brief Default constructor
*/
species::species()
{
	id = 0;
	sequence = "NRG";
	amount = 0;
	chargedMols = 0;
	concentration = 0;
	age = 0;
	reborns = 0;
	diffusionEnh = 1;
	soluble = SOLUBLE;
	complexCuttingPoint = 0;
	complexDegradationEnh = 0;
	evaluated = 1;
	catalyst_ID = 0;
	substrate_ID = 0;
	K_phospho = 0;
	energizable = 0;
	concentrationFixed = false;
	lastSpeciesEvaluated = 0;
	firstConcentration = 0;
 
}

/** \brief This constructor is used each time a new species is created (AMOUNT BASED)
 @param tmpID species identificator
 @param tmpSequence species sequence (e.g. ABABAABABA)
 @param tmpAmount species initial amount
 @param tmpDiffusionEnh Diffusion enhancement degree 
 @param tmpSoluble 1 if the species is soluble, 0 otherwise
 @param tmpComplexDegEnh complex dissociation kinetic constant
 @param tmpComplexCuttingPoint complex cutting point (catalyst-substrate)
 @param tmpEvalueted This parameter indicates whether the species has been already evalutad (i.e. all the catalysis of the species are instantiated)
 @param tmpVolume the volume is necessary to convert numbers in concentrations
 @param tmpK_phospho phosphorilation kinetic constant (in case of energy based simulations)
 @param tmpEnergizable this is a flag indicating whether or not the species is energizable
*/
species::species(acs_longInt tmpID, string tmpSequence, acs_longInt tmpAmount, acs_double tmpDiffusionEnh, 
				 acs_int tmpSoluble, acs_double tmpComplexDegEnh, acs_int tmpComplexCuttingPoint, 
				 acs_int tmpEvalueted, acs_double tmpVolume, acs_double tmpK_phospho, acs_int tmpEnergizable, acs_double tmpInflux_rate, acs_int tmpMaxLOut)
{
	id = tmpID;
	sequence = tmpSequence;
	amount = tmpAmount;
	chargedMols = 0;
	numToConc(tmpVolume);
	age = 0;
	reborns = 0;
	diffusionEnh =tmpDiffusionEnh;
	soluble = tmpSoluble;
	complexDegradationEnh = tmpComplexDegEnh;
	complexCuttingPoint = tmpComplexCuttingPoint;
	evaluated = tmpEvalueted;
	catalyst_ID = 0;
	substrate_ID = 0;
	K_phospho = tmpK_phospho;
	energizable = tmpEnergizable;
	
	// Set concentrationFixed propriety
	concentrationFixed = false;
	if((tmpInflux_rate == 0) && (tmpMaxLOut > 0))
	{
		if(sequence.length() <= tmpMaxLOut) concentrationFixed = true;
	}

    firstConcentration = 0; // 0 because the species is created during on run
    lastSpeciesEvaluated = 0;

}

/** \brief This constructor is used each time a new species is created (CONCENTRATION BASED)
 @param tmpID species identificator
 @param tmpSequence species sequence (e.g. ABABAABABA)
 @param tmpConcentration species initial concentration
 @param tmpDiffusionEnh Diffusion enhancement degree 
 @param tmpSoluble 1 if the species is soluble, 0 otherwise
 @param tmpComplexDegEnh complex dissociation kinetic constant
 @param tmpComplexCuttingPoint complex cutting point (catalyst-substrate)
 @param tmpEvalueted This parameter indicates whether the species has been already evalutad (i.e. all the catalysis of the species are instantiated)
 @param tmpVolume the volume is necessary to convert concentrations in numbers
 @param tmpK_phospho phosphorilation kinetic constant (in case of energy based simulations)
 @param tmpEnergizable this is a flag indicating whether or not the species is energizable
 */
species::species(acs_longInt tmpID, string tmpSequence, acs_double tmpConcentration, acs_double tmpDiffusionEnh, 
				 acs_int tmpSoluble, acs_double tmpComplexDegEnh, acs_int tmpComplexCuttingPoint, 
				 acs_int tmpEvalueted, acs_double tmpVolume, acs_double tmpK_phospho, acs_int tmpEnergizable, acs_double tmpInflux_rate, acs_int tmpMaxLOut)
{
	id = tmpID;
	sequence = tmpSequence;
	concentration = tmpConcentration;
	concToNum(tmpVolume);
	numToConc(tmpVolume);
	chargedMols = 0;
	age = 0;
	reborns = 0;
	diffusionEnh =tmpDiffusionEnh;
	soluble = tmpSoluble;
	complexDegradationEnh = tmpComplexDegEnh;
	complexCuttingPoint = tmpComplexCuttingPoint;
	evaluated = tmpEvalueted;
	catalyst_ID = 0;
	substrate_ID = 0;
	K_phospho = tmpK_phospho;
	energizable = tmpEnergizable;
	
	// Set concentrationFixed propriety
	concentrationFixed = false;
	if((tmpInflux_rate == 0) && (tmpMaxLOut > 0))
	{
		if(sequence.length() <= tmpMaxLOut) concentrationFixed = true;
	}

    firstConcentration = 0; // 0 because the species is created during on run
    lastSpeciesEvaluated = 0;
}

/** \brief This constructor is used when a new species is uploaded from file (TOTAL AMOUNT BASED)
 @param tmpID species identificator
 @param tmpSequence species sequence (e.g. ABABAABABA)
 @param tmpAmount species initial amount of molecules
 @param tmpReactions_constant ????
 */
species::species(acs_longInt tmpID, string tmpSequence, acs_longInt tmpAmount, acs_double tmpDiffusionEnh, 
				 acs_int tmpSoluble, acs_double tmpComplexDegEnh, acs_int tmpComplexCuttingPoint, 
				 acs_int tmpEvalueted, acs_double tmpAge, acs_int tmpReborns, acs_double tmpVolume, 
                 acs_longInt tmpNotUsedCatID, acs_longInt tmpNotUsedSubID, acs_double tmpK_phospho,
                 acs_int tmpEnergizable, acs_double tmpInflux_rate, acs_int tmpMaxLOut)
{

	id = tmpID;
	sequence = tmpSequence;
	amount = tmpAmount;
	chargedMols = 0;
    numToConc(tmpVolume); // The concetration is set within this function
	age = tmpAge;
	reborns = tmpReborns;
	diffusionEnh =tmpDiffusionEnh;
	soluble = tmpSoluble;
	complexDegradationEnh = tmpComplexDegEnh;
	complexCuttingPoint = tmpComplexCuttingPoint;
	evaluated = tmpEvalueted;
	catalyst_ID = tmpNotUsedCatID;
	substrate_ID = tmpNotUsedSubID;
	K_phospho = tmpK_phospho;
	energizable = tmpEnergizable;
	
	// Set concentrationFixed propriety
	concentrationFixed = false;
	if((tmpInflux_rate == 0) && (tmpMaxLOut > 0))
	{
		if(sequence.length() <= tmpMaxLOut) concentrationFixed = true;
	}

    firstConcentration = concentration; // If the species is loaded from file that's the very initial concentration
    lastSpeciesEvaluated = 0;
	
}

/** \brief This constructor is used when a new species is uploaded from file (CONCENTRATION BASED)
 @param tmpID species identificator
 @param tmpSequence species sequence (e.g. ABABAABABA)
 @param tmpAmount species initial amount of molecules
 @param tmpReactions_constant ????
 */
species::species(acs_longInt tmpID, string tmpSequence, acs_double tmpConcentration, acs_double tmpDiffusionEnh, 
				 acs_int tmpSoluble, acs_double tmpComplexDegEnh, acs_int tmpComplexCuttingPoint, 
				 acs_int tmpEvalueted, acs_double tmpAge, acs_int tmpReborns, acs_double tmpVolume, 
                 acs_longInt tmpNotUsedCatID, acs_longInt tmpNotUsedSubID, acs_double tmpK_phospho,
                 acs_double tmpKLoadConc, acs_int tmpEnergizable, acs_double tmpInflux_rate, acs_int tmpMaxLOut)

{

	id = tmpID;
	sequence = tmpSequence;
	concentration = tmpConcentration;
	concToNum(tmpVolume);
	numToConc(tmpVolume);
	chargedMols = round(tmpKLoadConc*tmpVolume*AVO);
	age = tmpAge;
	reborns = tmpReborns;
	diffusionEnh =tmpDiffusionEnh;
	soluble = tmpSoluble;
	complexDegradationEnh = tmpComplexDegEnh;
	complexCuttingPoint = tmpComplexCuttingPoint;
	evaluated = tmpEvalueted;
	catalyst_ID = tmpNotUsedCatID;
	substrate_ID = tmpNotUsedSubID;
	K_phospho = tmpK_phospho;
	energizable = tmpEnergizable;
	
	// Set concentrationFixed propriety
	concentrationFixed = false;
	if((tmpInflux_rate == 0) && (tmpMaxLOut > 0))
	{
		if(sequence.length() <= tmpMaxLOut) concentrationFixed = true;
	}

    firstConcentration = tmpConcentration; // If the species is loaded from file that's the very initial concentration
    lastSpeciesEvaluated = 0;
}

/** \brief This constructor is used when a species is randomly created (!!! NOT USED NOW)
 @version 0.1 (8 parameters)
 @param tmpID species identificator
 @param tmpSequence species sequence (e.g. ABABAABABA)
 @param tmpAmount species initial amount of molecules
 @param acs_double tmpDiffusionEnh Diffusione enhancement parameter
 @param acs_double tmpPrecipitationEnh Precipitation Enhancement parameters
 @param acs_double tmpComplexProb Probability to be a complex
 @param acs_double tmpMaxComplexDegKinetic max complex degradation constant
 @param MTRand& tmp_RandomGenerator random generator
 */
species::species(acs_longInt tmpID, string tmpSequence, acs_longInt tmpAmount, acs_double tmpDiffusionEnh, 
				 acs_int tmpSoluble, acs_double tmpComplexProb, acs_double tmpMaxComplexDegKinetic, 
				 MTRand& tmp_RandomGenerator, acs_double tmpVolume, acs_double tmpK_phospho, acs_int tmpEnergizable)
{
	cout << "four" << endl;

	id = tmpID;
	sequence = tmpSequence;
	amount = tmpAmount;
	numToConc(tmpVolume);
	chargedMols = 0;
	age = 0;
	reborns = 0;
	diffusionEnh =tmpDiffusionEnh;
	soluble = tmpSoluble;
	// RANDOM SELECTION BETWEEN SPECIES AND COMPLEX ACCORDING TO THE COMPLEX CREATION PROBABILITY
	acs_double dice = tmp_RandomGenerator();
	if(dice <= tmpComplexProb) //COMPLEX
	{
		complexDegradationEnh = getDoubleRandom(0.0, tmpMaxComplexDegKinetic, tmp_RandomGenerator);
		complexCuttingPoint = getIntRandom(2, sequence.length(), tmp_RandomGenerator) - 1;
		
	}else{ //SPECIES
		complexDegradationEnh = 0;
		complexCuttingPoint = 0;
	}
	evaluated = 0;
	K_phospho = tmpK_phospho;
	energizable = tmpEnergizable;
	concentrationFixed = false;
    firstConcentration = 0; // 0 because the species is created during on run
    lastSpeciesEvaluated = 0;
}

/** \brief This constructor is used to create a molecular complex
 @version 0.1 (10 paramters)
 @param tmpID species identificator
 @param tmpSequence species sequence (e.g. ABABAABABA)
 @param acs_double tmpDiffusionEnh Diffusione enhancement parameter
 @param acs_double tmpPrecipitationEnh Precipitation Enhancement parameters
 @param acs_double tmpComplexProb Probability to be a complex
 @param acs_double tmpMaxComplexDegKinetic max complex degradation constant
 @param MTRand& tmp_RandomGenerator random generator
 @param acs_int tmpCatalyst_ID Catalyst ID
 @param acs_int tmpSubstrate_ID substrate ID
 */
species::species(acs_longInt tmpID, string tmpSequence, acs_double tmpDiffusionEnh, 
				 acs_int tmpSoluble, acs_double tmpMaxComplexDegKinetic, 
				 acs_int tmpCuttingPoint, MTRand& tmp_RandomGenerator, acs_longInt tmpCatalyst_ID, 
				 acs_longInt tmpSubstrate_ID, acs_double tmpVolume, acs_double tmpK_phospho, acs_int tmpEnergizable)
{
	id = tmpID;
	sequence = tmpSequence;
	amount = 1;
	chargedMols = 0;
	numToConc(tmpVolume);
	age = 0;
	reborns = 0;
	diffusionEnh =tmpDiffusionEnh;
	soluble = tmpSoluble;
	complexDegradationEnh = getDoubleRandom(tmpMaxComplexDegKinetic, tmpMaxComplexDegKinetic, tmp_RandomGenerator);
	complexCuttingPoint = tmpCuttingPoint;
	evaluated = 0;
	catalyst_ID = tmpCatalyst_ID;
	substrate_ID = tmpSubstrate_ID;
	K_phospho = tmpK_phospho;
	energizable = tmpEnergizable;
	concentrationFixed = false;
    firstConcentration = 0; // 0 because the species is created during on run
    lastSpeciesEvaluated = 0;
}

/** Function to increment the total number of molecules belonging to this species
 */
void species::increment(acs_double tmpVolume)
{
    if(!concentrationFixed)
   {
       amount++;
       numToConc(tmpVolume);
   }
}

/** Function to decrement the total number of molecules belonging to this species
 */
void species::decrement(acs_double tmpVolume)
{
   if(!concentrationFixed)
   {
       amount--;
       numToConc(tmpVolume);
   }
}

/** \Funtion to charge a specific number of molecules
 */
bool species::setSpecificChargeMols(acs_int tmpMolsToCharge)
{
	bool chargeOK = false;
	chargedMols += tmpMolsToCharge;
	if (chargedMols > amount) {
		chargedMols = amount;
		chargeOK = true;
	}
	return chargeOK;
}

/** \Funtion to charge a specific number of molecules
 */
bool species::setChargeMols(acs_int tmpMolsToCharge)
{
	bool chargeOK = false;
	chargedMols = tmpMolsToCharge;
	if (chargedMols > amount) {
		chargedMols = amount;
		chargeOK = true;
	}
	return chargeOK;
}

/** \Funtion to charge molecules
*/
bool species::chargeMol()
{
	bool chargeOK = false;
	if (chargedMols < amount) {
		chargedMols++;
		chargeOK = true;
	}else{
		cout << "+: " << chargedMols << " tot: " << amount << endl;
	}

	return chargeOK;
}

/** \Funtion to uncharge molecules
 */
bool species::unchargeMol()
{
	bool unchargeOK = false;
	if (chargedMols > 0) {
		chargedMols--;
		unchargeOK = true;
	}
	return unchargeOK;
}

/** \Function return the concentration of the loaded molecules 
*/
acs_double species::getLoadedConcentration(acs_double tmpVolume)
{
	acs_double tmpConc = chargedMols / (AVO * tmpVolume);
	//cout << "Amount: " << amount << " - Conc:" << concentration << " - Loaded M:" << chargedMols << " -  Conc C:" << tmpConc << endl;
	return tmpConc;
}
/**
 \Function to insert the second substrate - and k reaction - to the species (complexes list)
 */
void species::insertSecSub(acs_longInt tmpID, acs_double tmpK, acs_longInt tmpCat)
{
    if(secondSubstrates.size() > 0)
    {
    	bool tmpFound = false;
    	for(acs_int z = 0; z < (acs_int)secondSubstrates.size(); z++)
    	{
    		if(secondSubstrates.at(z) == tmpID){tmpFound = true; break;}
    	}
    	if(!tmpFound){secondSubstrates.push_back(tmpID); kCond.push_back(tmpK); catalysisIfCpx.push_back(tmpCat);}
    }else{
    	secondSubstrates.push_back(tmpID);
    	kCond.push_back(tmpK);
    	catalysisIfCpx.push_back(tmpCat);
    }
}

/**
 * Function to show the species gillespie engagement
 */
void species::showGillEngagement()
{
	cout << "|- SPECIES " << id << " is engaged in Gillespie reactions: ";
	for(acs_int i = 0; i < gillespieEngagement.size(); i++)
	{
		cout << gillespieEngagement.at(i) << "\t";
	}
	cout << endl;
}
