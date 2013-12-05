/*
 *  reactions.cpp
 *  acsm2s
 *
 *  Created by Alessandro Filisetti on 31/03/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include "reactions.h"

/** \brief Constructor 
 @param tmpID reaction identificator
 @param tmpType condensation or cleavage
 @param tmpM_I product (if condensation) or substrates (if cleavage)
 @param tmpM_II product (if cleavage) or substrates (if condensation)
 @param tmpM_III product (if cleavage) or substrates (if condensation)
 @param tmpKeq equilibrium constant
 */
reactions::reactions(acs_longInt tmpID, acs_int tmpType, acs_longInt tmpM_I, acs_longInt tmpM_II, 
                     acs_longInt tmpM_III, acs_int tmpEvents, acs_int tmpEnergyType, acs_double tmpKspont)
{
	id_reaction = tmpID;
	type = tmpType;
	species_I = tmpM_I;
	species_II = tmpM_II;
	species_III = tmpM_III;
	events = tmpEvents;
	energyType = tmpEnergyType;	
	k_spont = tmpKspont;
}

