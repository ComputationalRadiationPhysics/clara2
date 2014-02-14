//***************************************************************************//
//  mainClara2.cc
//***************************************************************************//
//
//  Copyright: (C) 2014 by Alexander Debus, Richard Pausch
//  Email    : a.debus@hzdr.de
//  Authors  : Alexander Debus, Richard Pausch
//
//***************************************************************************//

///\file mainClara2.hh This file is included from both, GPT and Clara2. Thus, it has to be compatible to C *and* C++.

#pragma once

///\todo Hopefully temporally...
#ifndef NUMPARS_SUPPORTED
#define NUMPARS_SUPPORTED 10000000
#endif //NUMPARS_SUPPORTED

#ifndef OUTPUT_INTERVAL
#define OUTPUT_INTERVAL 100
#endif //OUTPUT_INTERVAL

#ifdef __cplusplus
extern  "C" {
#endif // __cplusplus

///The structure in which old information is kept
typedef struct {
  double GBx, GBy, GBz;
  double G;
  char stored;
} old_par;

extern int step;	//Index of the step
extern int pars_clara2;
extern int alive;  //Number of particles calculated using Clara2 (output)

///This variable keeps a second particles-array
extern old_par* pars_old;

#ifdef __cplusplus
} //This brace belongs to "extern C".
#endif // __cplusplus

//The cpp header follows
#ifdef __cplusplus

//------------------------------------------------------------------------//
//---- C++ header here ---------------------------------------------------//
//------------------------------------------------------------------------//
//------------------------------------------------------------------------//
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <cassert>
#include <cstdlib>
#include <fstream>

// Tell mainClara2 about the GPT data structures
#include <math.h>
#ifndef INCLUDED_ELEM
#define INCLUDED_ELEM
#include "elem.h"
#endif //INCLUDED_ELEM

using namespace std;

#include "vector.hpp"
#include "detector_dft.hpp"
#include "physics_units.hpp"
#include <omp.h>	// OpenMP

///The routines for Clara2 calculation.
class mainClara2 {

	public:
		///Constructor
		mainClara2();
		///Shutdown Clara2.
		~mainClara2();
		
		/* ------------ constants ------------------------------- */
		static const double omega_max = 1.2e21;   				/* maximum of ploted frequency Hz */
		static const double theta_max = 2.5/300;			  	/* maximum of theta in rad */
		static const unsigned int N_spectrum = 2048;			/* number of frequencies "omega"*/
		static const unsigned int N_theta = 120; 				/* number of directions in first angle "theta" */
		static const unsigned int N_phi = 2;       		 		/* number of directions in second angle "phi" */
		static const unsigned int N_direction = N_theta*N_phi; 	/* number of all directions */
		
		///Init Clara2. Has to be called before any other function.
		bool init();
		struct spectrum_container
		{
			std::complex<double> spectrum_x[N_spectrum];
			std::complex<double> spectrum_y[N_spectrum];
			std::complex<double> spectrum_z[N_spectrum];
		};
		static spectrum_container* all_spec;					/* All resulting spectral radiation amplitudes are saved here. */
	
    bool calculate(
					int numpar,							/*!< Number of entries in the gpt particle array. */
					double time,						/*!< The current GPT simulation time.*/
					double dt							/*!< The current GPT time step.*/
                  );
    
    void outputData();
    
    ///Returns a pointer to the global instance of mainClara2. If it does not exist, this function creates a new one.
    static mainClara2* getInstance() {return (instance == NULL) ? (instance = new mainClara2()) : instance;};
	
	private:
		static mainClara2* instance;
		int numtasks, rank;

		double omega[N_spectrum];								/* Array of frequencies to be calculated. */
		double theta[N_theta];									/* Array over theta angles */
		double phi[N_phi];										/* Array over phi angles */
		
};


//Mark the next functions to be compiled compatible to C
extern "C" {
//End of C++ part, the following is parsed by both
#endif // __cplusplus

//------------------------------------------------------------------------//
//---- Wrapper functions here --------------------------------------------//
//------------------------------------------------------------------------//
//------------------------------------------------------------------------//

	int clara2_init();

	int clara2_calculate(int numpar, double time, double dt);

	void outputData();

#ifdef __cplusplus
} //This brace belongs to "extern C".
#endif // __cplusplus

