/**
 * Copyright 2014-2016 Richard Pausch, Alexander Koehler
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


#include <iostream>
#include <cstdio>
#include <fstream>
#include <sstream>
#include "gzip_lib.hpp"
#include "settings.hpp"
#include "setFilename.hpp"


template<typename T, unsigned int N>
struct spectrum
{
  spectrum(void)
  {
    for(unsigned int i=0; i<N; ++i)
      spec[i] = (T) 0.0;
  }

  T spec[N];
};




int main(int argc, char * const argv[])
{

  using namespace std;
  const unsigned int N_direction = param::N_theta * param::N_phi;
  const unsigned int N_split = 8;

  spectrum<double, param::N_omega>* data = new spectrum<double, param::N_omega>[N_direction];


  // ----- run through all input files -------
  for(unsigned int index_files = param::index_files_first; index_files <= param::index_files_last; ++index_files)
  {
    char filename[param::N_char_filename];
    setFilename(filename, param::outputFileTemplate, index_files, param::N_char_filename);

    // ------ read input file ------
    if(param::ascii_output)
    {
      // ---- files loaded have ASCII  format -------
      FILE* pFile = fopen(filename, "r");
      if(pFile == NULL)
      {
        std::cout << "abort ascii file " << index_files
                  << "/" << param::index_files_last << std::endl;
        continue;
      }
      else
      {
        std::cout << "load ascii file " << index_files
                  << "/" << param::index_files_last << std::endl;
      }

      for(unsigned int j=0; j<N_direction; ++j)
      {
        double data_dump[N_split];
        for(unsigned int i=0; i< param::N_omega; i+=N_split)
        {
          if(fscanf(pFile, "%lf \t%lf \t%lf \t%lf \t%lf \t%lf \t%lf \t%lf",
                    &(data_dump[0]), &(data_dump[1]), &(data_dump[2]), &(data_dump[3]),
                    &(data_dump[4]), &(data_dump[5]), &(data_dump[6]), &(data_dump[7]))
             == (int)N_split)
          {
            for(unsigned int a =0; a<N_split; ++a)
              data[j].spec[i+a] += data_dump[a];
          }
          else
            throw "data";
        }
      }
      fclose(pFile);
    }
    else
    {
      // ---- binary input ---------
      const unsigned int N_double = param::N_omega * N_direction;
      double* data_input = new double[N_double];

      if(read_data(data_input, N_double*sizeof(double), filename) == 1)
      {
        std::cout << "abort binary file " << index_files
                  << "/" << param::index_files_last << std::endl;
        continue;
      }
      else
      {
        std::cout << "load binary file " << index_files
                  << "/" << param::index_files_last << std::endl;
      }

      for(unsigned int i=0; i<N_double; ++i)
        data[i/param::N_omega].spec[i%param::N_omega] += data_input[i];

      delete[] data_input;
    }

  }


  std::cout << "store data" << std::endl;
  for(unsigned int index_phi = 0; index_phi < param::N_phi; ++index_phi)
  {
    char output_filename[param::N_char_filename];
    setFilename(output_filename, param::output_pattern, index_phi, param::N_char_filename);

    std::ofstream output(output_filename);
    if(output.is_open())
    {
      for(unsigned j=index_phi*param::N_theta; j< (param::N_theta*(index_phi+1)); ++j)
      {
        for(unsigned i=0; i<param::N_omega; ++i)
        {
          output << data[j].spec[i] << " \t";
        }
        output << std::endl;
      }
      output.close();
    }
    else
    {
      std::cerr << "error: could not write output-file --> phi-index:" << index_phi << std::endl;
      throw "error output";
    }
  }

  delete[] data;


  std::cout << "deleting files" << std::endl;
  for(unsigned int index_files = param::index_files_first;
      index_files <= param::index_files_last;
      ++index_files)
  {
    char filename[param::N_char_filename];
    setFilename(filename, param::outputFileTemplate, index_files, param::N_char_filename);

    if(remove(filename) == 0)
      std::cout << "removed: " << filename << std::endl;
    else
      std::cerr << "error removing: " << filename << std::endl;
  }



  return 0;
}
