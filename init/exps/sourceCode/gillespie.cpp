/*
 *  gillespie.cpp
 *  acsm2s
 *
 *  Created by Alessandro Filisetti on 22/04/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */
 
 //prova2 svn

#include "gillespie.h"


/**
 Gillespie Obj constructor
 @version 1.0
 @param tmpRndDoubleGen randomGenerator reference
*/

gillespie::gillespie(acs_longInt tmpIDU, acs_int tmpIdReactionType, acs_double tmpScore, acs_longInt tmpMolI, 
					 acs_longInt tmpMolII, acs_longInt tmpMolIII, acs_longInt tmpMolIV, acs_longInt tmpIdReaction, 
					 acs_longInt tmpIdCatalysis)
{
	id_u = tmpIDU;
	idReactionType = tmpIdReactionType;
	score = tmpScore;
	molI = tmpMolI;
	molII = tmpMolII;
	molIII = tmpMolIII;	
	molIV = tmpMolIV;
	idReaction = tmpIdReaction;
	idCatalysis = tmpIdCatalysis;
	nrgSide = NOTHINGLOAD;

}

gillespie::gillespie(acs_longInt tmpIDU, acs_int tmpIdReactionType, acs_double tmpScore, acs_longInt tmpMolI, 
					 acs_longInt tmpMolII, acs_longInt tmpMolIII, acs_longInt tmpMolIV, acs_int tmpNRGside, acs_longInt tmpIdReaction, 
					 acs_longInt tmpIdCatalysis)
{
	id_u = tmpIDU;
	idReactionType = tmpIdReactionType;
	score = tmpScore;
	molI = tmpMolI;
	molII = tmpMolII;
	molIII = tmpMolIII;	
	molIV = tmpMolIV;
	idReaction = tmpIdReaction;
	idCatalysis = tmpIdCatalysis;
	nrgSide = tmpNRGside;
	
}