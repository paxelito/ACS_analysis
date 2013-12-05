/*
 *  commonFunctions.cpp
 *  acsm2s
 *
 *  Created by Alessandro Filisetti on 27/03/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

//prova

#include "commonFunctions.h"

/// This funtion returns a random position in a probability weight array of N elements
acs_longInt returnSelectionIdFromAWeightProbVector(acs_double *tmpArray, MTRand& tmpRandomGenerator, acs_int tmpRow)
{
        // tmpArray is not a vector. Access is made by []
	acs_double pnt = tmpRandomGenerator();
	bool foundFlag = false;
	int position = 0;
	int idReturned = position;
	while (foundFlag == false) {
            if (tmpArray[position] > pnt)
            {
                    idReturned = position;
                    foundFlag = true;
            }
            position++;
	}
	return idReturned;
}

/**
 Return position of a randomly selected element from a vector containing cumulative values for each element
 @version 1.0
 @param vector<acs_double>& tmpVector
 @param acs_double MAX VALUE contained within the QList (being a cumulative list this is the last value)
 @param MTRand& tmpRandomGenerator
 */
acs_longInt returnSelectionIdFromAWeightProbVector(vector<acs_double>& tmpVector, acs_double tmpMaxValue, MTRand& tmpRandomGenerator, acs_int tmpRow)
{

    acs_double pnt = tmpRandomGenerator();
    bool foundFlag = false;
    acs_longInt position = 0;
    acs_longInt idReturned = position;
    while (foundFlag == false) {
        try{
            if (tmpVector.at(position) > (pnt*tmpMaxValue))
            {
                idReturned = position;
                foundFlag = true;
            }
        }catch(exception& e)
        {
            cout << "Error in returnSelectionIdFromAWeightProbVector tmpVector.at(position)...used in row " << tmpRow << endl;
            cout << "Vector size " << tmpVector.size() << " - position " << position << endl;
            cerr << "exception caught: " << e.what() << endl;
            ExitWithError("returnSelectionIdFromAWeightProbVector","exception error thrown");
        }
        position++;
    }
    return idReturned;
}

/**
 Return position of a randomly selected element from a normalized vector containing cumulative values for each element
 @version 1.0
 @param vector<acs_double>& tmpQList
 @param MTRand& tmpRandomGenerator
 */
acs_longInt returnSelectionIdFromAWeightProbVectorAlreadyNormalized(vector<acs_double>& tmpVector, MTRand& tmpRandomGenerator)
{
    acs_double pnt = tmpRandomGenerator();
    bool foundFlag = false;
    acs_longInt position = 0;
    acs_longInt idReturned = position;
    acs_double tmpCumulativeProbability = tmpVector.at(0);

    while (foundFlag == false) {

        if(position >= (unsigned)tmpVector.size())
        {
            cout << "! ! ! ! ! ! ! ! ! FATAL ERROR ! ! ! ! ! ! ! !" << endl;
            cout << "idReturned: " << idReturned << " Position " << position << " is grater than the maximum number "
            << tmpVector.size()  << " cum prob: " << tmpCumulativeProbability<< " RN: " << pnt << endl;
            cout << "\t Suggestion >> Check that the cumulative sum of the probabilities is equal to 1" << endl;
            ExitWithError("returnSelectionIdFromAWeightProbQListAlreadyNormalized", "Error in return selection from weight prob Qlist");
        }

        if (tmpCumulativeProbability >= pnt)
        {
            idReturned = position;
            foundFlag = true;
        }

        if (position == (unsigned)tmpVector.size()-1) {
            idReturned = position;
            foundFlag = true;
        }
        try{
            tmpCumulativeProbability += tmpVector.at(position);
        }catch(exception&e)
        {
            cout<<"returnSelectionIdFromAWeightProbVectorAlreadyNormalized tmpVector.at(position)..."<<endl;
            cout<<"Vector size "<<tmpVector.size()<<"- position"<<position<<endl;
            cerr<<"exceptioncaught:"<<e.what()<<endl;
            ExitWithError("returnSelectionIdFromAWeightProbVectorAlreadyNormalized","exceptionerrorthrown");
        }


        position++;

    }
    return idReturned;
}

/**
 Return position of a LONG INT randomly selected element from a normalized vector containing cumulative values for each element
 @version 1.0
 @param vector<acs_double>& tmpVector
 @param MTRand& tmpRandomGenerator
 */
acs_longInt returnUniformSelection_LONG_IdFromVector(vector<acs_longInt>& tmpVector, MTRand& tmpRandomGenerator)
{
    acs_longInt totalSize = tmpVector.size() - 1;
    acs_longInt pnt = getIntRandom(0, totalSize, tmpRandomGenerator);
    return pnt;
}

acs_longInt getIntRandom(acs_longInt tmpFromNum, acs_longInt tmpToNum, MTRand& tmpRandomGenerator)
{

	if(tmpToNum < tmpFromNum)
	{
		cout << "tmpToNum (" << tmpToNum << ") is < of tmpFromNum (" << tmpFromNum << ") ... Cannot generate a number with a negative delta" << endl;
		exit(EXIT_FAILURE);
	}
	
	acs_longInt Result = tmpRandomGenerator() * (tmpToNum + 1 - tmpFromNum) + tmpFromNum;
	// if the random outcome in Result variable is out of the boundaries passed by parameters
	// the function exit with failure
	if((Result < tmpFromNum) || (Result > tmpToNum))
	{
		cout << "Result: " << Result << "\tGenerate a number out of boundaries" << endl;
		exit(EXIT_FAILURE);
	}
	return(Result);
}

acs_double getDoubleRandom(acs_double tmpFromNum, acs_double tmpToNum, MTRand& tmpRandomGenerator)
{
	// if parameters passed to this function are negative
	// this function exit with failure

	if((tmpFromNum < 0) || (tmpToNum < 0))
	{
		cout << "tmpFromNum: " << tmpFromNum << "\ttmpToNum: " << tmpToNum << "\tCannot generate a number with negative numbers" << endl;
		exit(EXIT_FAILURE);
	}
	// if tmpToNum is minor of tmpFromNum this 
	// function exit with failure
	if(tmpToNum < tmpFromNum)
	{
		cout << "tmpToNum (" << tmpFromNum << ") is < of tmpFromNum (" << tmpToNum << ") ... Cannot generate a number with a negative delta" << endl;
		exit(EXIT_FAILURE);
	}
	
	MTRand randomGen;
	acs_double Result = tmpRandomGenerator() * (tmpToNum - tmpFromNum) + tmpFromNum;	// if the random outcome in Result variable is out of the boundaries passed by parameters
	// the function exit with failure
	if((Result < tmpFromNum) || (Result > tmpToNum))
	{
		cout << "Result: " << Result << "\tGenerate a number out of boundaries" << endl;
		exit(EXIT_FAILURE);
	}
	return(Result);
}

/**
 Function to return a number from a poisson random distribution 
 */
acs_longInt random_poisson(acs_double tmpLambda, MTRand& tmpRandomGenerator)
{
	acs_double p=exp(-tmpLambda);
	acs_double g=p;
	acs_double u=tmpRandomGenerator();
	acs_longInt k=0;
	if(g > 0)
	{
		while(u>g)
		{
			p*=(tmpLambda/(acs_double)(++k));
			g+=p;
		}
	}
	return k;
}

/**
 Function to return a number from a binomial  distribution 
 */
acs_longInt random_binomial(acs_longInt n, acs_double tmpP, MTRand& tmpRandomGenerator)
{
//	acs_double t=tmpP/(1-tmpP);
//	acs_double u=tmpRandomGenerator();
//	acs_double p0=pow(double(1-tmpP),double(n));
//	acs_double g=p0;
//	acs_longInt k=0;
//	while (u>g)
//    {
//		p0*=t*(n-k)/(k+1);
//		g+=p0;
//		k++;
//    };
//	return k;
	acs_longInt k = 0;
	for(acs_int i = 0; i < n; i++)
	{
		if(tmpRandomGenerator() < tmpP)
			k++;
	}
	return k;
}

/**
 Function to round double numbers in integers 
*/
acs_double acsround(acs_double tmpX)
{
	// integer part
	acs_int integerPart = (acs_int)tmpX;
	acs_double decimalPart = tmpX - (acs_double)integerPart;
	if(decimalPart>=0.5){return ceil(decimalPart) + integerPart;}else{return floor(decimalPart) + integerPart;}	
}

/**
  Function to convert a decimal number in a binary string composed of 12 bit such
  Example -> input: 10, binary 1010 --> 000000001010
*/
string dec2bin(acs_int tmpInt){
    const int size=sizeof(tmpInt)*8;
    bool s=0;
    string res;
    for (int a=0;a<size;a++){
        bool bit=tmpInt>>(size-1);
        if (bit)
            s=1;
        if (s)
            res.push_back(bit+'0');
        tmpInt<<=1;
    }
    if (!res.size())
        res.push_back('0');

    // Add '0' char to complete the string of 12 characters
    while (res.size() < 12)
    {
          res = '0' + res;
    }

    return res;
}


/**
 Function to close the program after having en error
 */
bool ExitWithError(string strFunctionName, string strError) {
        cout << "Error in " << strFunctionName << endl;
        cout << "\t" << strError << endl;
        exit(EXIT_FAILURE);
        return true;
}

/**
  Function to split a string and save tokens in a vector
  */
vector<string> split(string str,const char* delim)
{
    char * writable = new char[str.size() + 1];
    std::copy(str.begin(), str.end(), writable);
    writable[str.size()] = '\0'; // don't forget the terminating 0

    char* token = strtok(writable,delim);

    vector<string> result;

    while(token != NULL)
    {
        result.push_back(token);
        token = strtok(NULL,delim);
    }

    delete[] writable;
    return result;
}

/**
 * Function to convert from string to boolean
 */
bool fromStrToBool(string const &string)
{
	return string[0] == '1';
}




