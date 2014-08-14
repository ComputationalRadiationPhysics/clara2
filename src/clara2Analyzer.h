/**
 * Copyright 2014 Alexander Debus, Lucas Clemente (flint)
 *
 * This file is part of Clara 2.
 *
 * Clara 2 is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Clara 2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Clara 2.
 * If not, see <http://www.gnu.org/licenses/>.
 */


#pragma once


#include <stdio.h>
#include <math.h>
#include <omp.h>



/* only include elem.h if not loaded previously */
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
