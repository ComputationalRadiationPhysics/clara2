/**
 * Copyright 2014-2016 Richard Pausch, Joy, Alexander Koehler
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
#include <omp.h>  // OpenMP

#include "single_direction.hpp"
#include "vector.hpp"
#include "physics_units.hpp"
#include "import_from_file.hpp"
#include "load_txt.hpp"
#include "include/input_output.hpp"
#include "settings.hpp"
#include "setFilename.hpp"


/**
 * function that calculates spectra in different directions for
 * a single particle trace
 *
 * @param trace_id a unique id which identifies the trajectory file
 *
 * @return error code
 **/
int all_directions(const unsigned int trace_id)
{

  using namespace std;


  /* ------ start time measurement ------------------ */
  struct timeval t1, t2;
  gettimeofday(&t1, NULL);


  /* ------------ constants ------------------------------- */
  const unsigned int N_direction = param::N_theta * param::N_phi; // number of all directions

  /* ---------- get trace ID ----------------- */

  /* check whether the id of the trace "trace_id" is larger than the given N_trace value */
  if(!(trace_id <= param::N_trace))
  {
    std::cout << "trace-ID is out of range (MAX = " << param::N_trace << ")" << std::endl;
    return 1;
  }


  /* ------- set up/compute all angles thetas ------------------ */
  double theta[param::N_theta];
  for(unsigned i=0; i< param::N_theta; ++i)
  {
    theta[i] = (double)i / param::N_theta * param::theta_max;
  }


  /* ------- set up/compute all angles phis ---------------------- */
  double phi[param::N_phi] = {0.0, 90.0};



  /* allocate memory for all spectra */
  struct spectrum_container
  {
    double spectrum[param::N_spectrum];
  };
  spectrum_container* all_spec = new spectrum_container[N_direction];


  /* compute the frequency array "omega" and fill spectra with zeros */
  double omega[param::N_spectrum];
  const double my_delta_omega = param::omega_max/param::N_spectrum;
  for(unsigned i=0; i<param::N_spectrum; ++i)
  {
    /* compute frequency */
    omega[i] = i * my_delta_omega;

    /* initialize spectra with zeros */
    for(unsigned j=0; j<N_direction; ++j)
      all_spec[j].spectrum[i] = 0.0;
  }


  /* ------- location of data ----------------------- */

  char filenameTrace[param::N_char_filename];
  setFilename(filenameTrace, param::traceFileTemplate, trace_id, param::N_char_filename);
  /* print out path name in order to check it in output files */
  std::cout << "check: filename: " << filenameTrace << std::endl;


  /* -------- load trace from file ------- */

  /* check if given file exists */
  if(!file_exists(filenameTrace))
    return 1;

  /* output to inform user that file is loaded */
  std::cout << "load file: " << filenameTrace << std::endl;

  /* create memory */
  const unsigned linenumber = linecounter(filenameTrace); /* get lines of data */
  one_line* data = new one_line[linenumber]; /* get memory for data */
  /* run function that fills data from file into "data": */
  load_txt(filenameTrace, linenumber, data);



  /* --------- calculate spectrum of one trace for all direction ---------- */


  /* in case of additional parallelization using OpenMP remove comment: */
  /* #pragma omp parallel for num_threads(4) schedule(dynamic, 1) */
  for(unsigned direction_index = 0; direction_index< N_direction; ++direction_index)
  {
    const double my_theta = theta[direction_index % param::N_theta];
    const double my_phi   = phi[direction_index/param::N_theta];
    printf("calculate direction: %4d -> theta: %3.5f , phi: %3.5f \n", direction_index, my_theta, my_phi);

    /*
     * compute the spectra for a single direction
     * and trow an error if something goes wrong
     */
    if((single_direction(data, linenumber, omega, all_spec[direction_index].spectrum, param::N_spectrum, my_theta, my_phi))!=0)
    {
      std::cerr << "error occurred in single_direction function" << std::endl;
      throw "error in single_direction function";
    }
  }



  /* ------- outputfile ------------------------------ */
  /* allocate memory for name of output file */
  char outputfilename[param::N_char_filename];
  setFilename(outputfilename, param::outputFileTemplate, trace_id, param::N_char_filename);
  /* print name of output file */
  std::cout << "check: output-filename: " << outputfilename << std::endl;


  /* --- file output ------------- */

  /* store spectral data either as binary or as ascii data */
  if(param::ascii_output)
  {
    /* ---- ASCII output file ------------------------ */
    ofstream my_output(outputfilename); /* create file */
    if(my_output.is_open()) /* check if it is open */
    {
      for(unsigned j=0; j<N_direction; ++j) /* for all directions */
      {
        for(unsigned i=0; i<param::N_spectrum; ++i) /* for all frequencies */
        {
          /* print spectral data separated by tabs */
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
  else /* in case binary output was chosen */
  {
    /* ----- binary output file -------- */
    /* allocate memory to store spectral and directional data at once */
    double* output_data = new double[param::N_spectrum*N_direction];

    /*
     * "output_index" is the counter for the array behind "output_data"
     * it is initialized with zero and increased by one for each stored value
     */
    for(unsigned j=0, output_index=0; j<N_direction; ++j) /* for all directions */
    {
      for(unsigned i=0; i<param::N_spectrum; ++i) /* for all frequencies */
      {
        /* fill output-data-container with computed spectral data */
        output_data[output_index] = all_spec[j].spectrum[i];

        /* increase "output_index" by one to address new memory
         *  location in next loop cycle */
        output_index++;
      }
    }
    /* store (uncompressed) collected data using clara2's input_output.hpp */
    store_data(output_data,
               param::N_spectrum*N_direction*sizeof(double),
               outputfilename);

    delete[] output_data; /* free allocated memory after output */
  }


  delete[] all_spec; /* free spectral data */
  delete[] data; /* free traces */



  /* ---- end time measurement and output ------------ */
  /* print used compute time for information */
  gettimeofday(&t2, NULL);
  const long runtime = (t2.tv_sec - t1.tv_sec)*1000000 + (t2.tv_usec - t1.tv_usec);
  std::cout << "time needed: "
            << (double)runtime * 1.0e-6
            << " sec" << std::endl;
  std::cout << "time per direction (and trace): "
            << (double)runtime / (double)(N_direction) * 1.0e-6
            << " sec" << std::endl;


  /* return NULL after all went well for this trajectory */
  return 0;
}
