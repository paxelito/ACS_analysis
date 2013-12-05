/*
 *  reactions.h
 *  acsm2s
 *
 *  Created by Alessandro Filisetti on 31/03/09.
 *  Modified by Alessandro Filisetti on 19/04/2011.
 *  Copyright 2011 Alessandro Filisetti. All rights reserved.
 *  Some comments here
 */

#ifndef REACTIONS
#define REACTIONS

#include "acs_headers.h"
#include "commonFunctions.h"

class reactions
{
private:
    acs_longInt id_reaction;
    acs_int type;               /**< 1 cleavage 0 condensation */
    acs_longInt species_I;      /**< species[species_II] + species[species_III] */
    acs_longInt species_II;
    acs_longInt species_III;
    acs_int events;             /**< Reaction counter */
    acs_int energyType;         /**< Energetic Boolean Function */
    acs_double k_spont;			/**< If spontanoues reactions are turned on this is the kinetic constant of the spontanous reaction */
	
	
public:
	reactions(acs_longInt tmpID, acs_int tmpType, acs_longInt tmpM_I, acs_longInt tmpM_II, 
                          acs_longInt tmpM_III, acs_int tmpEvents, acs_int tmpEnergyType, acs_double tmpKspont);
	~reactions() {}
	
	// GETTING FUNCTIONS
	acs_longInt getID()const{return id_reaction;}
	acs_int getType()const{return type;}
	acs_longInt getSpecies_I()const{return species_I;}
	acs_longInt getSpecies_II()const{return species_II;}
	acs_longInt getSpecies_III()const{return species_III;}
	acs_int getEvents()const{return events;}
	acs_int getEnergyType()const{return energyType;}
	acs_double getKspont()const{return k_spont;}

	// OPERATIVE FUNCTIONS
    void updateTotEvents(){events++;}
    void resetEventsCounter(){events = 0;}
};

#endif
