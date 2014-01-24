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



#include "all_directions.hpp"




#include <iostream>
#include <string>
#include <cassert>
#include <cstdlib>

#include <stdlib.h>
#include <sys/time.h>
#include <omp.h>	// OpenMP

#include "single_trace.hpp"

#include "vector.hpp"
#include "analytical_solution.hpp"
#include "physics_units.hpp"
#include "string_manipulation.hpp"


#include "import_from_file.hpp"
#include "load_txt.hpp"
#include "gzip_lib.hpp"





int all_directions(const unsigned int trace_id, const char arg[])
{

  using namespace std;


  // ------ start time measurement ------------------
  struct timeval t1, t2;
  gettimeofday(&t1, NULL);



  // ------------ constants -------------------------------
  const double omega_max = 3.0e19;  // maximum of ploted frequency Hz
  const double theta_max = 1.14594939; // maximum of theta in degree
  const unsigned int N_spectrum = 2048; // number of frequencies
  const unsigned int N_theta = 120;     // number of directions in theta
  const unsigned int N_phi = 2;         // number of directions in phi
  const unsigned int N_trace = 2000;    // maximum number of traces
  const unsigned int N_direction = N_theta*N_phi; // number of all directions

  // ---------- get trace ID -----------------

  if(!(trace_id <= N_trace))
    {
      std::cout << "trace-ID is out of range (MAX = " << N_trace << ")" << std::endl;
      return 1;
    }

  // -------- get store info -----------
  bool asci_output;
  std::string store_str = arg;
  if(!store_str.compare("asci"))
    {
      asci_output = true;
      std::cout << "ASCI output" << std::endl;
    }
  else if(!store_str.compare("binary"))
    {    
      asci_output = false;
      std::cout << "binary output" << std::endl;
    }
  else
    {
      std::cerr << "2nd argument needs to be binary or asci" << std::endl;
      throw "bin_asci";  
    }





  // ------- create all thetas --------------------------
  double theta[N_theta];
  for(unsigned i=0; i< N_theta; ++i)
    {
      theta[i] = (double)i / N_theta * theta_max;
    }


  // ------- create all phis ----------------------------
  double phi[N_phi] = {0.0, 90.0};




  /* ----- create spectrum for all ------ */
  double omega[N_spectrum];
  struct spectrum_container
  {
    double spectrum[N_spectrum];
  };
  spectrum_container* all_spec = new spectrum_container[N_theta*N_phi];
  

  const double my_delta_omega = omega_max/N_spectrum;
  for(unsigned i=0; i<N_spectrum; ++i)
    {
      omega[i] = i * my_delta_omega;
      for(unsigned j=0; j<N_direction; ++j)
	all_spec[j].spectrum[i] = 0.0;
    }



  // ------- location of data -------------------------
  const char directory[] = "/net/cns/projects/HPLsim/xray/debus/ELBEThomson/basicRun2/";
  const char prefix[] = "trace_";
  const char postfix[] = ".txt";
  char filename[256];
  if(sprintf(filename, 
	     "%s%s%04d%s",
	     directory,
	     prefix,
	     trace_id,
	     postfix) > 254)
    {
      std::cerr << "buffer  to small!!! " << std::endl;
      throw "Buffer to small!";
    }
  std::cout << "check: filename: " << filename << std::endl;


 
  // ------------------ load trace from file ---------------------------   
    
    if(!file_exists(filename)) 
      return 1;
    std::cout << "load file: " << filename << std::endl;
    const unsigned linenumber = linecounter(filename); // get lines of data    
    one_line* data = new one_line[linenumber]; // get memory for data
    load_txt(filename, linenumber, data);


  


  // --------- calculate spectrum of one trace for all direction ----------


    //  #pragma omp parallel for num_threads(4) schedule(dynamic, 1) 
  for(unsigned direction_index = 0; direction_index< N_direction; ++direction_index)
    {
      const double my_theta = theta[direction_index % N_theta];
      const double my_phi   = phi[direction_index/N_theta];
      printf("calculate direction: %4d -> theta: %3.5f , phi: %3.5f \n", direction_index, my_theta, my_phi);

      if((single_trace(data, linenumber, omega, all_spec[direction_index].spectrum, N_spectrum, my_theta, my_phi))!=0)
	{
	  std::cerr << "error occured in single_trace function" << std::endl;
	  throw "error in single_trace function";
	}
    }




  // ------- outputfile -------------------------------
  char outputfilename[256];
  if(sprintf(outputfilename, 
	     "my_spectrum_trace%06d.dat",
	     trace_id) > 254)
    {
      std::cerr << "buffer  to small!!! " << std::endl;
      throw "Buffer to small!";
    }
  std::cout << "check: output-filename: " << outputfilename << std::endl;



  // --- file output -------------

  if(asci_output)
    {
      // ---- ASCI output file ------------------------//
      ofstream my_output(outputfilename);
      if(my_output.is_open())
	{
	  for(unsigned j=0; j<N_direction; ++j)
	    {
	      for(unsigned i=0; i<N_spectrum; ++i)
		{
		  my_output << all_spec[j].spectrum[i] << " \t";
		}
	      my_output << std::endl;
	    }
	  my_output.close();
	}
      else
	{
	  std::cerr << "error writing output" << std::endl;
	  throw "error output";
	}
    }
  else
    {
      // ----- binary output file --------
      double* output_data = new double[N_spectrum*N_direction];
      for(unsigned j=0, output_index=0; j<N_direction; ++j)
	{
	  for(unsigned i=0; i<N_spectrum; ++i)
	    {
	      output_data[output_index] = all_spec[j].spectrum[i];
	      output_index++;
	    }
	}
      store_data(output_data, N_spectrum*N_direction*sizeof(double), 
		 outputfilename);

      delete[] output_data;
    }







  delete[] all_spec;
  delete[] data;


    

  // ---- end time measurement and output ------------
  gettimeofday(&t2, NULL);
  const long runtime = (t2.tv_sec - t1.tv_sec)*1000000 + (t2.tv_usec - t1.tv_usec);
  std::cout << "time nedded: " << (double)runtime * 1.0e-6 
	    << " sek" << std::endl;
  std::cout << "time per direction (and trace): " 
	    << (double)runtime / (double)(N_theta*N_phi) * 1.0e-6
	    << " sek" << std::endl;




  return 0;
}



