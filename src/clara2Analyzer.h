//***************************************************************************//
//  clara2Analyzer.h
//***************************************************************************//
//
//  Copyright: (C) 2014 by Alexander Debus
//  Email    : a.debus@hzdr.de
//  Authors  : Alexander Debus, Lucas Clemente (GPT Template)
//
//***************************************************************************//

#pragma once

#include <stdio.h>
#include <math.h>
#include <omp.h>

#ifndef INCLUDED_ELEM
#define INCLUDED_ELEM
#include "elem.h"
#endif //INCLUDED_ELEM

#include "mainClara2.hh"

///\todo Hopefully temporally...
#ifndef NUMPARS_SUPPORTED
#define NUMPARS_SUPPORTED 10000000
#endif //NUMPARS_SUPPORTED

#ifndef OUTPUT_INTERVAL
#define OUTPUT_INTERVAL 100
#endif //OUTPUT_INTERVAL

int step;	//Index of the step
int pars_clara2;
int alive;  //Number of particles calculated using Clara2 (output)

///This variable keeps a second particles-array
old_par* pars_old;
