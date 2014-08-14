/**
 * Copyright 2014 Alexander Debus, Richard Pausch
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



#include "mainClara2.hh"

//---------------------------------------------------------------------------
//--- Compability functions for GPT (C++ to C) ------------------------------
//---------------------------------------------------------------------------

extern "C" int clara2_init() {
	return mainClara2::getInstance()->init();
}

extern "C" int clara2_calculate(int numpar, double time, double dt)
{
	mainClara2::getInstance()->calculate(numpar, time, dt);
	return true;
}

extern "C" void outputData()
{
	mainClara2::getInstance()->outputData();
}


//---------------------------------------------------------------------------
//--- Real Clara2 class ---------------------------------------------------------
//---------------------------------------------------------------------------

mainClara2::spectrum_container* mainClara2::all_spec;

mainClara2::mainClara2()
{
	if (instance) printf("Clara2 already instantiated. The application could show undefined behavior!\n");
	instance = this;
}

mainClara2::~mainClara2()
{
	// Finalize
	delete[] all_spec; /* free spectral data */
	instance = 0;
	printf("Exit done. Clara2 is down...\n");
}

mainClara2* mainClara2::instance = 0;

//mainClara2::spectrum_container all_spec[mainClara2::N_direction];

bool mainClara2::calculate(int numpar, double time, double dt) {

	alive = 0;
	pars_clara2 = 0; //Number of particles calculated using Clara2 (output)

	printf("OpenMP 4.0: %d devices available.\n", omp_get_num_devices() );
	//Iterate over all particles
	///\todo Check whether gcc3 (GPT) is compatible to gcc4 (OpenMP)...
	//#pragma omp parallel for private(i)
	#pragma omp target data map( all_spec ) map( to: pars[0:numpar], pars_old[0:numpar]  )
	for (int i = 0; i < numpar; i++)
	{
		//printf("OpenMP: %d executes on host device.\n", omp_is_initial_device() );
		//Check if we have a problem
		if (pars[i].ID >= NUMPARS_SUPPORTED) {
			fprintf( stderr, "pars[%d].ID=%d >= NUMPARS_SUPPORTED=%d\n", i, pars[i].ID, NUMPARS_SUPPORTED);
		}
		//Is the particle still alive?
		if (!pars[i].alive) continue;
		alive++;

		//Is the particle in old_pars?
		if (pars_old[pars[i].ID].stored == 0)
		{
		  continue;
		}
		//If the particle has not moved, abort
		if ((pars[i].GBr[0] == pars_old[pars[i].ID].GBx) && (pars[i].GBr[1] == pars_old[pars[i].ID].GBy) && (pars[i].GBr[2] == pars_old[pars[i].ID].GBz))
		{
		  continue;
		}
		pars_clara2++;

		double* GBr_pre = (double*)malloc(3*sizeof(double));
		GBr_pre[0] = pars_old[pars[i].ID].GBx;
		GBr_pre[1] = pars_old[pars[i].ID].GBy;
		GBr_pre[2] = pars_old[pars[i].ID].GBz;
		double G_pre = pars_old[pars[i].ID].G;
		
		//// Needed for output in a particle's first step...
		//// This is probably outdated, because copyArray() insures, that pars[i].G is stored in pars_old.
		//if (pars_old[pars[i].ID].stored == 0) pars_old[pars[i].ID].G = pars[i].G;
		

		/* -------- Arrange data in memory for DFT------- */
		R_vec location = R_vec( pars[i].Wr[0] , pars[i].Wr[1] , pars[i].Wr[2] );
		R_vec beta = R_vec( pars[i].GBr[0]/pars[i].G , pars[i].GBr[1]/pars[i].G , pars[i].GBr[2]/pars[i].G );
		R_vec beta_dot = R_vec( 0.0 , 0.0 , 0.0 );
		if (dt > 0.0)
		{
			beta_dot[0]=(pars[i].GBr[0]/pars[i].G - GBr_pre[0]/G_pre) / dt;
			beta_dot[1]=(pars[i].GBr[1]/pars[i].G - GBr_pre[1]/G_pre) / dt;
			beta_dot[2]=(pars[i].GBr[2]/pars[i].G - GBr_pre[2]/G_pre) / dt;
		}

		/* --------- calculate spectral contribution of one timestep for all direction ---------- */

		/* in case of additional parallelsation using OpenMP uncomment this: */
		#pragma omp target
		#pragma omp teams distribute parallel for schedule(dynamic, 1)
		for(unsigned direction_index = 0; direction_index< N_direction; ++direction_index)
		{
			//printf("OpenMP 4.0: %d executes on host device.\n", omp_is_initial_device() );
			const double my_theta = theta[direction_index % N_theta];
			const double my_phi   = phi[direction_index/N_theta];
			
			/* Create, initialize DFT detector */
			R_vec looking_vector = R_vec(std::cos(my_theta) , std::sin(my_theta) * std::cos(my_phi) , std::sin(my_theta) * std::sin(my_phi) );
			Detector_dft dft(looking_vector, dt, N_spectrum, omega);
			
			//printf("calculate direction: %4d -> theta: %3.5f , phi: %3.5f \n", direction_index, my_theta, my_phi);
			
			/* DFT calculation of the spectral contribution of the current timestep. Then add the (complex-valued vector) result to the global spectral data all_spec */
			dft.add_to_spectrum(location, beta, beta_dot, time);
			for (unsigned i=0; i<N_spectrum; ++i)
			{
					all_spec[direction_index].spectrum_x[i]+=dft.spektrum[i][0];
					all_spec[direction_index].spectrum_y[i]+=dft.spektrum[i][1];
					all_spec[direction_index].spectrum_z[i]+=dft.spektrum[i][2];
			}
			
		}
		
		free(GBr_pre);
	}

	return true;
}

void mainClara2::outputData()
{
	/* ------- outputfile ------------------------------ */
	/* allocate memory for name of output file */
	char outputfilename[256];

	//Physical prefactor of the "Jackson formular"
	const double factor = util::square(phy::q)/(16.*util::cube(M_PI)*phy::epsilon_0*phy::c);
	
	// Auxiliary variables
	double intensity_x, intensity_y, intensity_z, intensity_total;
	
	/* fill output file for each trace_id based on template */
	if(sprintf(outputfilename, "my_spectrum.dat") > 254)
    {
      std::cerr << "buffer  to small!!! " << std::endl;
      throw "Buffer to small!";
    }
	/* print name of output file */
	std::cout << "check: output-filename: " << outputfilename << std::endl;
	
	/* ---- ASCI output file ------------------------ */
	ofstream my_output(outputfilename,ofstream::out); /* create file */
	if(my_output.is_open()) /* check if it is open */
	{
		for(unsigned j=0; j<N_direction; ++j) /* for all directions */
		{
			for(unsigned i=0; i<N_spectrum; ++i) /* for all frequencies */
			{
				/* Calculate the modulus square of the resulting complex amplitudes and include the physical prefactors */
				intensity_x = std::abs(all_spec[j].spectrum_x[i]); intensity_x *= intensity_x*factor;
				intensity_y = std::abs(all_spec[j].spectrum_y[i]); intensity_y *= intensity_y*factor;
				intensity_z = std::abs(all_spec[j].spectrum_z[i]); intensity_z *= intensity_z*factor;
				intensity_total = intensity_x + intensity_y + intensity_z ;
				/* print spectral data separated by tabs */
				my_output << intensity_total << " \t";
			}
			/* separate each direction by a newline */
			my_output << std::endl;
		}
		/* close output file */
		my_output.close();
	}
	else /* if output file is not open create a warning */
	{
		std::cerr << "error writing output" << std::endl;
		throw "error output";
	}
  
}

bool mainClara2::init(void)
{
	printf("Launched Clara2 Radiation-Analyzer.\n");

	/** the function 'start_array()' calls the parallelization procedure
	* defined in 'parallel_jobs.h'. The variables 'numtasks' and 'rank'
	* are set by this function and allow to use the parallel distrubuted 
	* indexes to be used.
	**/
	//if(start_array(&numtasks, &rank) != 0) return 1;

	/* ------- set up/compute all angles thetas ------------------ */
	for(unsigned i=0; i< N_theta; ++i)
	{
		theta[i] = (double)i / N_theta * theta_max;
	}

	/* ------- set up/compute all angles phis ---------------------- */
	phi[0] = 0.0;
	phi[1] = M_PI/2;

	/* allocate memory for all spectra */
	//mainClara2::spectrum_container all_spec[mainClara2::N_direction]
	all_spec = new spectrum_container[N_direction];
  
	/* compute the frequency array "omega" and fill spectra with zeros */
	const double my_delta_omega = omega_max/N_spectrum;
	for(unsigned i=0; i<N_spectrum; ++i)
	{
		/* compute frequency */
		omega[i] = i * my_delta_omega; 
		/* initialise spectra with zeros */
		for(unsigned j=0; j<N_direction; ++j)
		{
			all_spec[j].spectrum_x[i] = std::complex<double>(0.0,0.0);
			all_spec[j].spectrum_y[i] = std::complex<double>(0.0,0.0);
			all_spec[j].spectrum_z[i] = std::complex<double>(0.0,0.0);
		}
	}

	return true;
}
