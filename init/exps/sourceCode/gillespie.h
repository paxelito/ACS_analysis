/*! \class gillespie 
 * \author Alessandro Filisetti
 * \version 0.1
 * \date 2009-04-21
 */

#ifndef GILLESPIE
#define GILLESPIE

#include "acs_headers.h"

class gillespie
{
private:
	acs_longInt id_u; /**<  ID */
	acs_int idReactionType; /**<  Reaction type (go to acs_header.h for more informations) */
	acs_double score; /**<  Single Gillespie propensity */
	acs_longInt molI; /**<  According to the specific reaction type molI,II,III and IV represents different species (catalyst rather than substract...) */
	acs_longInt molII;
	acs_longInt molIII;
	acs_longInt molIV;
	acs_longInt idReaction; /**<  reactions object foreing key */
	acs_longInt idCatalysis; /**<  catalysis object foreing key */
	acs_int nrgSide; /**<  catalystLoaded, SubstrateLoaded, BothLoaded, NOenergy */
public: 
	gillespie();
	gillespie(acs_longInt tmpIDU, acs_int tmpIdReactionType, acs_double tmpScore, 
			  acs_longInt tmpMolI, acs_longInt tmpMolII, acs_longInt tmpMolIII, 
			  acs_longInt tmpMolIV, acs_longInt tmpIdReaction, acs_longInt tmpIdCatalysis);
	gillespie(acs_longInt tmpIDU, acs_int tmpIdReactionType, acs_double tmpScore, 
			  acs_longInt tmpMolI, acs_longInt tmpMolII, acs_longInt tmpMolIII, 
			  acs_longInt tmpMolIV, acs_int tmpNRGside, acs_longInt tmpIdReaction, acs_longInt tmpIdCatalysis);
	~gillespie(){}
	
	// GET FUNCTIONS
	acs_longInt getID()const{return id_u;}
	acs_int getIdReactionType()const{return idReactionType;}
	acs_double getScore()const{return score;}
	acs_longInt getMolI()const{return molI;}
	acs_longInt getMolII()const{return molII;}
	acs_longInt getMolIII()const{return molIII;}
	acs_longInt getMolIV()const{return molIV;}
	acs_longInt getReactionID()const{return idReaction;}
	acs_longInt getCatalysisID()const{return idCatalysis;}
	acs_int getNRGside()const{return nrgSide;}
	


};


#endif
