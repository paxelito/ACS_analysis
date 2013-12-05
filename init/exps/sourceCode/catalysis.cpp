/*
 *  catalysis.cpp
 *  acsm2s
 *
 *  Created by Alessandro Filisetti on 16/04/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

//prova
#include "catalysis.h"


/**	catalysis class constructor (FROM FILE)
 * 
 * \version 0.1
 * \date 2010-03-16
 */
catalysis::catalysis(acs_longInt tmpCatId, acs_longInt tmpCat, acs_longInt tmpRctId, acs_longInt tmpAmount, 
					 acs_double tmpKass, acs_double tmpKdiss, acs_double tmpK_cpx, acs_int tmpCpxTarget)
{
	if((tmpCpxTarget < 1) || (tmpCpxTarget > 2))
	{
		cout << "ID: " << tmpCatId
			 << " || Catalyst: " << tmpCat
			 << " || reaction_ID: " << tmpRctId
			 << " || totalAmount: " << tmpAmount
			 << " || Kdiss: " << tmpKdiss
			 << " || Kass: " << tmpKass
			 << " || k_cpx: " << tmpK_cpx
			 << " || cpxTarget: " << tmpCpxTarget << endl;
		ExitWithError("catalysis::catalysis","Wrong Complex Target");
	}
	ID = tmpCatId;				// CATALYSIS ID
	catalyst = tmpCat;			// CATALYST ID
	reaction_ID = tmpRctId;		// REACTION ID
	totalAmount = tmpAmount;	// TOTAL AMOUNT OF THESE CATALYSIS 
	Kdiss = tmpKdiss;			// CLEAVAGE KINETIC CONSTANT
	Kass = tmpKass;				// CONDENSATION KINETIC CONSTANT
	k_cpx = tmpK_cpx;			// COMPLEX FORMATION KINETIC CONSTANT
	cpxTarget = tmpCpxTarget;	// Catalyst Complex species Target
	//TR k_cpxDiss = tmpK_cpxDiss;	// COMPLEX DISSOCIATION KINETIC CONSTANT
}
