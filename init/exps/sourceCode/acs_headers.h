/*
 *  acs_headers.h
 *  acsSIM
 *
 *  Created by Alessandro Filisetti on 09/03/09.
 *  Copyright 2009 __European Centre for Living Technology__. All rights reserved.
 *
 */

//prova2 cin.ignore().get();

#ifndef ACS_HEADERS
#define ACS_HEADERS

// LIBRARIES

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <sstream>
#include <vector>
#include <deque>
#include <cmath>
#include <cstdlib>
#include <clocale> // Probably this library can be removed after QT removing process will be finished
#include <ctime>
#include "mtrand.h"


// TYPEDEF

typedef long double	acs_double;
typedef unsigned long int acs_longInt;
typedef unsigned int acs_int;

//----------
// CONSTANTS
//----------

#define __SOFTVERSION__ "4.8b20131125.61"

// DEBUG LEVELS

#define MINIMAL_PROMPT -1
#define RUNNING_VERSION 0
#define SMALL_DEBUG 1
#define MEDIUM_DEBUG 2
#define HIGH_DEBUG 3
#define FINDERRORDURINGRUNTIME -10
#define COMPLEXSTUFF -20
#define GILLESPIESTUFF -30

// RANDOM INITIALIZATION
#define RANDOMRANGE random()

// INITIAL SPECIES CONCENTRATION DISTRIBUTION
#define PROPORTIONALMOLECULEAMOUNT 1
#define UNIFORMMOLECULEAMOUNT 2
#define INVPROPORTIONALMOLECULEAMOUNT 3

// POSSIBLE REACTIONS
#define CONDENSATION 0
#define CLEAVAGE 1
#define COMPLEXFORMATION 2
#define COMPLEXDEGRADATION 3
#define SPECIESDECAY 4
#define PHOSPHORILATION 5
#define ENDO_CLEAVAGE 6
#define ENDO_CONDENSATION 7
#define ENDO_COMPLEXFORMATION 8
#define ENERGYEFFLUX 9
#define SPONTANEOUS_CONDENSATION 10
#define SPONTANEOUS_CLEAVAGE 11

// SOLUBILITY
#define SOLUBLE 1
#define PRECIPITATED 0

// ENERGY
#define ESOERGONIC 1
#define ENDOERGONIC 0
#define CLEAVAGEBASED 1
#define CONDENSATIONBASED 0
#define ENERGYBASED 1
#define ENERGYFREE 0
#define TRUENRG '1'
#define FALSENRG '0'
#define ENERGIZABLE 1
#define NOTENERGIZABLE 0

#define SUBSTRATELOAD 0
#define CATALYSTLOAD 1
#define BOTHLOAD 2
#define COMPLEXLOAD 3
#define NOTHINGLOAD 4

// REACTION CREATION TYPE
#define NEWREACTIONS 1
#define UPGRADEREACTIONS 0

// REATIONS MISC
#define NOSPONTANEOUS 0

// SPECIAL NUMBERS
#define NEP 2.7182818284590452353602874
#define AVO 6.02214179e+23

// MINIMAL MULTIPLICATOR FOR MINIMAL REACTION TIME
#define MINIMALRCTTIMEMULTI 100

using namespace std;

#endif
