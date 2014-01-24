/**
 * Copyright 2014 Richard Pausch
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




#include "single_trace.hpp"

#include "interpolation.hpp"
#include "discrete.hpp"
#include "run_through_data.hpp"


int single_trace(const one_line* data, const unsigned int linenumber,
		 const double* all_omega, double* all_spectrum, const unsigned N_all_spec,
		 const double theta_offset, const double phi_offset)
{
  ////////////////////////////////////////////////////////////
  // CHANGE THIS: 


  //const char target[] = "../Thomson_scattering/thomson_trajectory_gamma=02_a0=2.dat";
  const unsigned N_angle_theta = 1;
  const unsigned N_angle_phi   = 1;


  //////////////////////////////////////////////////////////////


  
    
    /* ---------- time steps : delta_t ---------------------- */

    
    double stepwidth = (data[6].intern_data[6] -data[5].intern_data[6])/* *1.0E-15 */;
    //std::cout << "delta_t : " << stepwidth << std::endl;
    //std::cout << "lines   : " << linenumber << std::endl << std::endl;

    // the equi-distant time step aka: delta_t 
    




    /* ----------------------- detectors ------------------------ */


    // creating different angles:
    R_vec looking_vector[N_angle_theta * N_angle_phi];
    double angle_theta[N_angle_theta];
    double angle_phi[N_angle_phi];

    for(unsigned i=0; i< N_angle_theta; ++i)
      {
	angle_theta[i] = (double(i)/double(N_angle_theta)) * M_PI * 0.5 + theta_offset/180.0 * M_PI;
      }

    for(unsigned i=0; i< N_angle_phi; ++i)
      {
	angle_phi[i] = (double(i)/double(N_angle_phi))*1.0 * M_PI + phi_offset/180.0 * M_PI;
      }

    for (unsigned a=0; a<N_angle_theta; ++a)
      {
	for (unsigned b=0; b<N_angle_phi; ++b)
	  {
	    looking_vector[a*N_angle_phi+b] 
	      = R_vec(std::cos(angle_theta[a]) ,
		      std::sin(angle_theta[a])*std::cos(angle_phi[b])  , 
		      std::sin(angle_theta[a])* std::sin(angle_phi[b]) );
	  }
      }

    //std::cout << "signle-trace: theta= " << angle_theta[0] 
    //	      << " phi= " << angle_phi[0]
    //	      << " vector= " << looking_vector[0] << std::endl;


    //std::cout << "Done: angles" << std::endl;




    /* -------- FFT ------------- */


    Detector_fft* detector_fft[N_angle_theta*N_angle_phi];

    // create FFT detectors
    for(unsigned i=0; i<N_angle_theta*N_angle_phi; ++i)
      {
      detector_fft[i] = new Detector_fft((looking_vector[i]), linenumber );
      }
    //std::cout << "created FFT Detectors\n" << flush;


    // analyse trajectory
    run_through_data(data, linenumber, N_angle_theta*N_angle_phi, detector_fft); 

    // calculate spectra FFT 
    for(unsigned k=0; k<N_angle_theta*N_angle_phi; ++k)
      {
	//	std::cout << " calculate FFT spektrum: " << k+1 << "/" << N_angle_theta*N_angle_phi << std::endl;
	(*detector_fft[k]).calc_spectrum();
      }


#if 0

    /* ------------------ DFT  ------------------- */

    Detector_dft* detector_dft[N_angle_theta * N_angle_phi];

    // create DFT detectors (using omega[i] from FFT)
    for(unsigned i=0; i<N_angle_theta*N_angle_phi ; ++i)
      {    
	unsigned dummy_N =  (*detector_fft[i]).half_frequency();
	//double dummy_omega_max = (*detector_fft[i]).get_spectrum(dummy_N-1, 0);
	double* dummy_omega = (*detector_fft[i]).frequency;
	detector_dft[i] = new Detector_dft((looking_vector[i]), stepwidth, 
					   dummy_N, dummy_omega /*_max*/);
      }
    //std::cout << "created DFT Detectors\n" << flush;

    // analyse trajectory
    run_through_data(data, linenumber, N_angle_theta*N_angle_phi, detector_dft); 

    // calculate spectra DFT
    for(unsigned k=0; k< N_angle_theta*N_angle_phi ; ++k)
      {
	//	std::cout << " calculate DFT spektrum: " << k+1 << "/" << N_angle_theta*N_angle_phi << std::endl;
	(*detector_dft[k]).calc_spectrum();
      }

#endif


    /* ------ combine spectra ---------------- */

    //    std::cout << "  --> single ANGLE: theta: " << angle_theta[0] << "   phi: " << angle_phi[0] << std::endl;

    interpolation_int(detector_fft[0], all_omega, all_spectrum, N_all_spec);




    /* ------- spectrum oputput to file ------ */
#if 0

    if(angle_phi[0] == 0.0 && (angle_theta[0] == 0.0 || angle_theta[0] == 0.5 * 2.98622)) // last point compare double is critical
      {
	char buffer[10];
	sprintf(buffer,"%3.5f", angle_theta[0]);
	

	//std::cout << "store spectrum: " << theta_offset << " ?= " << angle_theta[0] ;

	string dummy = buffer; //name;
	//ofstream energy_dat_dft("energy_dft.dat");
	//ofstream energy_dat_fft(("energy_fft" + dummy + ".dat").c_str());
	//ofstream spectrum_dft("spectrum_dft.dat");
	ofstream spectrum_fft(("spectrum_fft" + dummy + ".dat").c_str());
	unsigned location;

	for(unsigned a=0; a < N_angle_theta; ++a)
	  {
	    for (unsigned b=0; b< N_angle_phi; ++b)
	      {
		location = a*N_angle_phi + b;
		//energy_dat_fft << (*detector_fft[location]).energy() << '\t';
		//energy_dat_dft << (*detector_dft[location]).energy() << '\t';
		for (unsigned i=0; i< (*detector_fft[0]).half_frequency(); ++i)
		  {
		    //spectrum_dft << angle_theta[a] << "\t";
		    //spectrum_dft << angle_phi[b] << "\t";
		    //spectrum_dft << (*detector_dft[location]).get_spectrum(i, 0) << "\t";
		    //spectrum_dft << (*detector_dft[location]).get_spectrum(i, 1) << std::endl;
		
		    spectrum_fft << angle_theta[a] << "\t";
		    spectrum_fft << angle_phi[b] << "\t";
		    spectrum_fft << (*detector_fft[location]).get_spectrum(i, 0) << "\t";
		    spectrum_fft << (*detector_fft[location]).get_spectrum(i, 1) << std::endl;
		
		  }
	      }
	    //energy_dat_fft << std::endl;
	    //energy_dat_dft << std::endl;
	  }
      }

#endif
       
    //std::cout << "Done \n \n------------------------\n";


    for(unsigned i=0; i<N_angle_theta*N_angle_phi; ++i)
      {
	//delete detector_dft[i];
	delete detector_fft[i];
      }

    return 0;
}


bool file_exists(const char *filename)
{
  ifstream ifile(filename);
  return ifile;
}
