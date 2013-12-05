/** \class catalysis 
 * \brief CATALYSIS class
 *
 * This class contains catalysis proprieties and methods
 * \author Alessandro Filisetti
 * \version 0.1
 * \date 2009-04-16
 */

#ifndef CATALYSIS
#define CATALYSIS

#include "acs_headers.h"
#include "commonFunctions.h"

class catalysis
{
private:
        acs_longInt ID;         /**<  catalysis ID */
        acs_longInt catalyst;   /**<  catalyst ID */
        acs_longInt reaction_ID;/**<  reaction ID, this reaction is catalysed by catalyst ID */
        acs_longInt totalAmount;/**<  counter of reaction occurences */
	
	//Kinetic parameters 
        acs_double Kass;	/**<  Associative constant (CONDENSATION) */
        acs_double Kdiss;	/**<  dissociative constant (CLEAVAGE) */
        acs_double k_cpx;	/**<  Kinetic constant for the complex formation */
        //TR acs_double k_cpxDiss;	/**<  Kinetic constant for the complex dissociation */
        acs_int cpxTarget; /**<  species forming the complex with the catalyst in reaction structure, mol_II (1) or mol_III (2) */
	
public:
	catalysis();
	catalysis(acs_longInt tmpCatId, acs_longInt tmpCat, acs_longInt tmpRctId, acs_longInt tmpAmount, 
			  acs_double tmpKass, acs_double tmpKdiss, acs_double tmpK_cpx, acs_int tmpCpxTarget);
	~catalysis(){};
	
	//GETTING FUNCTIONS
        acs_longInt getCatId()const{return ID;} // Primary Key
        acs_longInt getCat()const{return catalyst;} // Catalyst ID (Foreign Key of species object)
        acs_longInt getReactionID()const{return reaction_ID;} // Reaction ID (Foreign Key of reaction object)
        acs_longInt getTotAmount()const{return totalAmount;} // Catalysis counter
        acs_double getKass()const{return Kass;} // Association kinetic constant A_C + B -> AB + C
        acs_double getKdiss()const{return Kdiss;} // dissociation Kinetic Constant A_C -> A + C (A_C is a temporary complex)
        acs_double getK_cpx()const{return k_cpx;} // complex formation kinetic constant A + C -> A_C
        //acs_double getK_cpxDiss()const{return k_cpxDiss;} // This parameter in now specific for the single species
        acs_int getCpxTarget()const{return cpxTarget;} // complex target.
	
	// OPERATIVE FUNCTION
        void updateTotAmount(){totalAmount++;} // Increase the counter of this catalysis
        void resetEventsCounter(){totalAmount = 0;}
	
};

#endif
