/** \class commonFunctions 
 *  \brief This class contains all the common function of the system
 *
 *  This class contains all the functions useful in general
 *  \authors alessandro filisetti
 *  \date 2011/12/10
 *  \version 1.0
 */

#ifndef COMMONFUNCTIONS
#define COMMONFUNCTIONS

#include "acs_headers.h"

//FUNCTIONS
acs_longInt returnSelectionIdFromAWeightProbVector(acs_double *tmpArray, MTRand& tmpRandomGenerator, acs_int tmpRow = 0);
acs_longInt returnSelectionIdFromAWeightProbVector(vector<acs_double>& tmpVector, acs_double tmpMaxValue, MTRand& tmpRandomGenerator, acs_int tmpRow = 0);
acs_longInt returnSelectionIdFromAWeightProbVectorAlreadyNormalized(vector<acs_double>& tmpVector, MTRand& tmpRandomGenerator);
acs_longInt returnUniformSelection_LONG_IdFromVector(vector<acs_longInt>& tmpVector, MTRand& tmpRandomGenerator);

acs_longInt getIntRandom(acs_longInt tmpFromNum, acs_longInt tmpToNum, MTRand& tmpRandomGenerator);
acs_double getDoubleRandom(acs_double tmpFromNum, acs_double tmpToNum, MTRand& tmpRandomGenerator);
acs_longInt random_poisson(acs_double tmpLambda, MTRand& tmpRandomGenerator);
acs_longInt random_binomial(acs_longInt n, acs_double tmpP, MTRand& tmpRandomGenerator);
acs_double acsround(acs_double tmpX);
string dec2bin(acs_int tmpInt);
vector<string> split(string str,const char* delim);
bool ExitWithError(string strFunctionName, string strError);
bool fromStrToBool(string const &string);

#endif
