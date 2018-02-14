/**
 * Copyright 2014-2018 Richard Pausch, Alexander Koehler
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


#include <fstream>
#include "include/input_output.hpp"
#include "settings.hpp"
#include "setFilename.hpp"


/** data structure to hold spectral data
  * that is automatically initialized with zeros
  */
template<typename T, unsigned int N>
struct spectrum
{
  /** constructor
    * initialize data with zeros
    */
  spectrum(void)
  {
    /* initialize all spectral data with zeros */
    for(unsigned int i=0; i<N; ++i)
      spec[i] = (T) 0.0;
  }

  T spec[N]; /* data */
};



/** process all spectra and generate total (incoherent) spectra
  *
  * Reads all spectra from each trace and combines them incoherently
  * to a total spectra for each phi value.
  *
  * @return int error code: 0 -> successful
  */
int main()
{
  using namespace std;

  /* compute total number of observation direction */
  const unsigned int N_direction = param::N_theta * param::N_phi;

  /* spectral container to hold total spectrum (for each direction) */
  spectrum<double, param::N_omega>* data = new spectrum<double, param::N_omega>[N_direction];


  /* -----------------------------
     READ SPECTRA FROM INPUT FILES
     -----------------------------*/
  for(unsigned int index_files = param::index_files_first; index_files <= param::index_files_last; ++index_files)
  {
    /* set file name for spectra from trace with id index_files */
    char filename[param::N_char_filename];
    setFilename(filename, param::outputFileTemplate, index_files, param::N_char_filename);

    /* ----- read input file ------ */
    if(param::ascii_output) /* if files loaded have ASCII format */
    {
      FILE* pFile = fopen(filename, "r"); /* file handler */

      if(pFile == NULL) /* file handler doe not point to file - abort this file */
      {
        std::cout << "abort ascii file " << index_files
                  << "/" << param::index_files_last << std::endl;
        continue;
      }
      else /* file handler points to valid file */
      {
        std::cout << "load ascii file " << index_files
                  << "/" << param::index_files_last << std::endl;
      }

      /* process all directions */
      for(unsigned int j=0; j<N_direction; ++j)
      {
        /* chunk size (in doubles) for data loading for ascii files */
        const unsigned int N_split = 8;
        double data_dump[N_split]; /* temporal data container */

        /* go through all frequencies */
        for(unsigned int i=0; i< param::N_omega; i+=N_split)
        {
          /* read data from ascii file and convert to double */
          if(fscanf(pFile, "%lf \t%lf \t%lf \t%lf \t%lf \t%lf \t%lf \t%lf",
                    &(data_dump[0]), &(data_dump[1]), &(data_dump[2]), &(data_dump[3]),
                    &(data_dump[4]), &(data_dump[5]), &(data_dump[6]), &(data_dump[7]))
             == (int)N_split)
          {
            for(unsigned int a =0; a<N_split; ++a)
              data[j].spec[i+a] += data_dump[a]; /* add spectral data to container */
          }
          else
            throw "data"; /* in case data reading fails - throw error */
        }
      }
      fclose(pFile);  /* close file handler */
    }
    else /* if files loaded have binary format */
    {
      /* create double array for entire spectra (of one file) */
      const unsigned int N_double = param::N_omega * N_direction;
      double* data_input = new double[N_double];

      /* read binary data at once */
      if(read_data(data_input, N_double*sizeof(double), filename) == 1)
      {
        /* if read failed - abort this file */
        std::cout << "abort binary file " << index_files
                  << "/" << param::index_files_last << std::endl;
        continue;
      }
      else /* if read was successful - print verbose output */
      {
        std::cout << "load binary file " << index_files
                  << "/" << param::index_files_last << std::endl;
      }

      /* add entire spectrum to all spectrum container
         (for all directions and frequencies) */
      for(unsigned int i=0; i<N_double; ++i)
        data[i/param::N_omega].spec[i%param::N_omega] += data_input[i];

      delete[] data_input; /* delete allocated memory */
    }

  }
  /* ---------------------------------
     END READ SPECTRA FROM INPUT FILES
     ---------------------------------*/


  /* ----------------------------------
     STORE INCOHERENT RADIATION SPECTRA
     ----------------------------------*/

  std::cout << "store data" << std::endl; /* verbose output */

  /* store spectra over theta  for each phi */
  for(unsigned int index_phi = 0; index_phi < param::N_phi; ++index_phi)
  {
    /* set output file name */
    char output_filename[param::N_char_filename];
    setFilename(output_filename, param::output_pattern, index_phi, param::N_char_filename);

    std::ofstream output(output_filename); /* file handler */
    if(output.is_open()) /* file successfully created */
    {
      /* for each theta do: */
      for(unsigned j=index_phi*param::N_theta; j< (param::N_theta*(index_phi+1)); ++j)
      {
        /* for each omega do: */
        for(unsigned i=0; i<param::N_omega; ++i)
        {
          output << data[j].spec[i] << " \t";
        }
        output << std::endl;
      }
      output.close();
    }
    else /* file creation failed - throw error */
    {
      std::cerr << "error: could not write output-file --> phi-index:" << index_phi << std::endl;
      throw "error output";
    }
  }

  delete[] data;

  /* --------------------------------------
     END STORE INCOHERENT RADIATION SPECTRA
     --------------------------------------*/


  /* -------------------------------------
     REMOVE SPECTRA FROM INDIVIDUAL TRACES
     -------------------------------------*/

  std::cout << "deleting files" << std::endl; /* verbose output */

  /* for each file previously loaded: */
  for(unsigned int index_files = param::index_files_first;
      index_files <= param::index_files_last;
      ++index_files)
  {
    /* set file name */
    char filename[param::N_char_filename];
    setFilename(filename, param::outputFileTemplate, index_files, param::N_char_filename);

    /* try deleting and give verbose output about success or failure */
    if(remove(filename) == 0)
      std::cout << "removed: " << filename << std::endl;
    else
      std::cerr << "error removing: " << filename << std::endl;
  }

  /* -----------------------------------------
     END REMOVE SPECTRA FROM INDIVIDUAL TRACES
     -----------------------------------------*/

  return 0;
}
