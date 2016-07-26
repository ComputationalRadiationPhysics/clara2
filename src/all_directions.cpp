/**
 * Copyright 2014-2016 Richard Pausch
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




/**
 * function that calculates spectra in different directions for
 * a single particle trace
 *
 * @param trace_id a unique id which which the trajctopry file 
 *                 can be identified 
 * @param arg  a string telling wether a "binary" or "ascii" 
 *             output should be used
 * @return error code
 **/
int all_directions(const unsigned int trace_id, const char arg[])
{

  using namespace std;


  /* ------ start time measurement ------------------ */
  struct timeval t1, t2;
  gettimeofday(&t1, NULL);



  /* ------------ constants ------------------------------- */
  const double omega_max = 3.0e19;      /* maximum of ploted frequency Hz */
  const double theta_max = 1.14594939;  /* maximum of theta in degree */
  const unsigned int N_spectrum = 2048; /* number of frequencies "omega"*/
  const unsigned int N_theta = 120;     /* number of directions in first angle "theta" */
  const unsigned int N_phi = 2;         /* number of directions in second angle "phi" */
  const unsigned int N_trace = 2000;    /* maximum number of traces */
  const unsigned int N_direction = N_theta*N_phi; // number of all directions

  /* ---------- get trace ID ----------------- */

  /* check whether the id of the trace "trace_id" is larger than the given N_trace value */
  if(!(trace_id <= N_trace))
    {
      std::cout << "trace-ID is out of range (MAX = " << N_trace << ")" << std::endl;
      return 1;
    }

  /* -------- get store info ----------- */
  /* since output can be stored as binary file or as ascii file, here the selected option
   * is checked or an error is thrown in case the selction was wrong */
  bool ascii_output;
  std::string store_str = arg;
  if(!store_str.compare("ascii"))
    {
      ascii_output = true;
      std::cout << "ASCII output" << std::endl;
    }
  else if(!store_str.compare("binary"))
    {    
      ascii_output = false;
      std::cout << "binary output" << std::endl;
    }
  else
    {
      std::cerr << "2nd argument needs to be binary or ascii" << std::endl;
      throw "bin_ascii";  
    }





  /* ------- set up/compute all angles thetas ------------------ */
  double theta[N_theta];
  for(unsigned i=0; i< N_theta; ++i)
    {
      theta[i] = (double)i / N_theta * theta_max;
    }


  /* ------- set up/compute all angles phis ---------------------- */
  double phi[N_phi] = {0.0, 90.0};




  /* allocate memory for all spectra */
  struct spectrum_container
  {
    double spectrum[N_spectrum];
  };
  spectrum_container* all_spec = new spectrum_container[N_theta*N_phi];
  


  /* compute the frequency array "omega" and fill spectra with zeros */
  double omega[N_spectrum];
  const double my_delta_omega = omega_max/N_spectrum;
  for(unsigned i=0; i<N_spectrum; ++i)
    {
      /* compute frqueny */
      omega[i] = i * my_delta_omega; 

      /* initialise spectra with zeros */
      for(unsigned j=0; j<N_direction; ++j)
        all_spec[j].spectrum[i] = 0.0; 
    }



  /* ------- location of data ----------------------- */

  /* TO DO: SIMPLIFY THIS BY USING SPRINTF() */
  /* set directory where to find the data: */
  const char directory[] = "/net/cns/projects/HPLsim/xray/debus/ELBEThomson/basicRun2/";
  /* set name of trajectory file before index appears in file name: */
  const char prefix[] = "trace_";
  /* set name of file after index is used */
  const char postfix[] = ".txt";

  /* join all parts together and store path to file in "filename" */
  char filename[256];
  if(sprintf(filename, 
	     "%s%s%04d%s",
	     directory,
	     prefix,
	     trace_id,
	     postfix) > 254)
    {
      /* throw warning when buffer is to small for path name */
      std::cerr << "buffer  to small!!! " << std::endl;
      throw "Buffer to small!";
    }
  /* print out path naem in order to check it in output files */
  std::cout << "check: filename: " << filename << std::endl;


 
  /* -------- load trace from file ------- */
    
  /* check if given file exists */
  if(!file_exists(filename)) 
    return 1;

  /* output to inform user that file is loaded */
  std::cout << "load file: " << filename << std::endl;

  /* create memory */
  const unsigned linenumber = linecounter(filename); /* get lines of data */
  one_line* data = new one_line[linenumber]; /* get memory for data */
  /* run function that fills data from file into "data": */
  load_txt(filename, linenumber, data); 


  


  /* --------- calculate spectrum of one trace for all direction ---------- */


  /* in case of additional parallelsation using OpenMP uncomment this: */
  /* #pragma omp parallel for num_threads(4) schedule(dynamic, 1) */
  for(unsigned direction_index = 0; direction_index< N_direction; ++direction_index)
    {
      const double my_theta = theta[direction_index % N_theta];
      const double my_phi   = phi[direction_index/N_theta];
      printf("calculate direction: %4d -> theta: %3.5f , phi: %3.5f \n", direction_index, my_theta, my_phi);

      /* 
       * compute the spectra for a single direction 
       * and trow an error if something goes wrong
       */
      if((single_trace(data, linenumber, omega, all_spec[direction_index].spectrum, N_spectrum, my_theta, my_phi))!=0)
        {
          std::cerr << "error occured in single_trace function" << std::endl;
          throw "error in single_trace function";
        }
    }




  /* ------- outputfile ------------------------------ */
  /* allocate memory for name of output file */
  char outputfilename[256];

  /* fill output file for each trace_id based on template */
  if(sprintf(outputfilename, 
	     "my_spectrum_trace%06d.dat",
	     trace_id) > 254)
    {
      std::cerr << "buffer  to small!!! " << std::endl;
      throw "Buffer to small!";
    }
  /* print name of output file */
  std::cout << "check: output-filename: " << outputfilename << std::endl;



  /* --- file output ------------- */

  /* store spectral data either as binary or as ascii data */
  if(ascii_output)
    {
      /* ---- ASCII output file ------------------------ */
      ofstream my_output(outputfilename); /* create file */
      if(my_output.is_open()) /* check if it is open */
        {
          for(unsigned j=0; j<N_direction; ++j) /* for all directions */
            {
              for(unsigned i=0; i<N_spectrum; ++i) /* for all frequencies */
                {
                  /* print spectral data seperated by tabs */
                  my_output << all_spec[j].spectrum[i] << " \t";
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
  else /* in case binary output was choosen */
    {
      /* ----- binary output file -------- */
      /* allocate memory to store spectral and directional data at once */
      double* output_data = new double[N_spectrum*N_direction];

      /* 
       * "output_index" is the counter for the array behind "output_data" 
       * it is inizialsed with zero and increased by one for each stored value
       */
      for(unsigned j=0, output_index=0; j<N_direction; ++j) /* for all directions */
        {
          for(unsigned i=0; i<N_spectrum; ++i) /* fora ll frequencies */
            {
              /* fill output-data-container with computed spectral data */ 
              output_data[output_index] = all_spec[j].spectrum[i];

              /* increase "output_ibdex" by one to address new memory 
               *  location in next loop cycle */
              output_index++; 
            }
        }
      /* store (uncompressed) collected data using clara2's gzib_lib.hpp */
      store_data(output_data, N_spectrum*N_direction*sizeof(double), 
                 outputfilename);
      
      delete[] output_data; /* free allocated memory after output */
    }




  delete[] all_spec; /* free spectral adat */
  delete[] data; /* free traces */


    

  /* ---- end time measurement and output ------------ */
  /* print used compute time for information */
  gettimeofday(&t2, NULL);
  const long runtime = (t2.tv_sec - t1.tv_sec)*1000000 + (t2.tv_usec - t1.tv_usec);
  std::cout << "time nedded: " << (double)runtime * 1.0e-6 
	    << " sek" << std::endl;
  std::cout << "time per direction (and trace): " 
	    << (double)runtime / (double)(N_theta*N_phi) * 1.0e-6
	    << " sek" << std::endl;


  /* return NULL after all went well for this trajectory */
  return 0;
}



