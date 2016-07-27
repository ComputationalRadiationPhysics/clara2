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


/**
 * This is a no longer used program's code that
 * that was needed when a single MPI task calculated the
 * radiation for only a single direction but for several 
 * particle traces. It collected all the results for 
 * different directions and created a single matrix-like
 * output that contained the radiation spectra from all
 * trajectories and all directions.
 *
 * SHOULD THIS STILL BE INCLUDED IN THE CODE? - ISSUE #17
 **/




#include <iostream>
#include <string>
#include <cstdlib>
#include <fstream>
#include <cmath>


/** TO DO: this should be in a separate file - ISSUE #15 **/

/**
 * This function checks if a file exists on the hard drive.
 *
 * @param filename string containing the path and filename to be checked
 * @return Returns true if file exists, otherwise false.
 **/
bool file_exists(const char *filename)
{
  std::ifstream infile(filename);
  return infile;
}



/**
 * This program combines files with spectra for different directions 
 * to a single files containing all spectra for directions.
 *
 * @return only 0 is return - no error code yet  
 **/
int main()
{

  /* verbose output to inform user that the program is running */
  std::cout << "Convert single file data to matrix: " << std::flush;

  /* ----------- set parameters ------------ */
  /* in this section, all parameters are set that need to be changed 
   * in order to set up the collection process. */

  const unsigned N_omega = 2048; /* number of frequencies */
  /* To DO: could be read from linenumber input file or from a 
   *  general param file  associated with the simulation - ISSUE #8 */

  double omega[N_omega]; /* allocate memory for reading omega values in */



  /* compute directions: theta and phi */
  /* The real values of theta and phi are needed because they are used
   * for naming the files. */



  /* compute directions: theta and phi */
  /* The real values of theta and phi are needed because they are used
   * for naming the files. */

  const unsigned N_theta = 120; /* number of different theta angles */
  /* To DO: this can not be read from input file but from a 
   *  general param file  associated with the simulation - ISSUE #8 */

  double theta[N_theta]; /* allocate memory for directions */

  /* compute directions (theta): */
  double theta_max =  1.14594939; /* = 1.0/gamma in degree */
  for(unsigned i=0; i<N_theta; ++i)
    {
      theta[i] = (double)i / N_theta * theta_max;
    }

  const unsigned N_phi = 2; /* number of different phi angles */
  /* To DO: this can not be read from input file but from a 
   *  general param file  associated with the simulation - ISSUE #8 */

  /* set phi angles */
  double phi[N_phi] = {0.0, 90.0};


  /* allocate memory for collecting spectra for different directions */
  double data_phi_0[N_theta][N_omega]; /* matrix data  for phi=0 degrees */
  double data_phi_90[N_theta][N_omega]; /* matrix data  for phi=90 degrees */



  /* -------- end set parameters -------- */




  /* -------- begin load data -------- */

  
  char buffer[10]; /* allocate memory for direction-id */

  /* run through all directions */
  for(unsigned i=0; i< N_phi*N_theta; ++i)
    {
      /* string with correct format for theta: */
      sprintf(buffer, "%3.5f", theta[i%N_theta]);
      std::string theta_st = buffer; 

      /* string with correct format for phi: */
      sprintf(buffer, "%3.5f", phi[i/N_theta]);
      std::string phi_st = buffer; 

      /* create file name for single direction (theta and phi) */
      std::string filename = "my_output_theta=" + theta_st + "_phi=" + phi_st + ".dat";

      /* check if input file exists and load it*/
      if(!(file_exists(filename.c_str())))
        std::cout << "  ERROR: file " << filename << " does not exists" << std::endl;
      else /* file exists */
        {
          /* connect to file */
          std::ifstream data(filename.c_str()); 

          /* verbose information to users which file is loaded */ 
          std::cout << "loading file: " << filename << std::endl;

          /* branch for different files: */
          if(i==0) /* with the fist file omega values are collected */
            {
              /* read omega AND spectra data from file */
              for(unsigned j=0; j<N_omega; ++j)
                data >> omega[j] >> data_phi_0[i][j];
            }
          else if(i<N_theta) /* for all other files of phi=0 */
            {
              double omega_temp; /* variable to temporarily store omega */ 
              for(unsigned j=0; j<N_omega; ++j)
                {
                  /* read spectra and omega (again) */
                  data >> omega_temp >> data_phi_0[i][j];

                  /* check if previously stored omegas are equal to 
                   * those stored in the other files. If they are not
                   * equal write warning to screen. (but do not stop 
                   * program) */
                  if(omega_temp != omega[j])
                    std::cout << "  ERROR: frequency discrepancy: " << filename
                              << " - " << omega[j] << " != " << omega_temp << std::endl;
                }
            }
          else /* for all files with phi=90 --> same procedure */
            {
              double omega_temp; /* temporary variable for omega */
              /* go through data in files (see above's code */
              for(unsigned j=0; j<N_omega; ++j)
                {
                  data >> omega_temp >> data_phi_90[i%N_theta][j];
                  if(omega_temp != omega[j])
                    std::cout << "  ERROR: frequency discrepancy: " << filename
                              << " - " << omega[j] << " != " << omega_temp << std::endl;
                }
            }
          data.close(); /* close the file of the current direction */
        }
    } 
  

  /* -------- end load data --------- */


  /* store matrix-like spectral data for first phi 
   * tabs separate the values of different frequencies
   * and newlines separate different directions. */
  std::ofstream output_0("matrix_phi_0.dat");
  for(unsigned i=0; i<N_theta; ++i) /* all theta */
    {
      for (unsigned j=0; j<N_omega; ++j) /* all omega */
        output_0 << data_phi_0[i][j] << " \t";
      output_0 << std::endl;
    }


  /* store matrix-like spectral data for second phi 
   * tabs separate the values of different frequencies
   * and newlines separate different directions. */
  std::ofstream output_90("matrix_phi_90.dat");
  for(unsigned i=0; i<N_theta; ++i) /* all theta */
    {
      for (unsigned j=0; j<N_omega; ++j) /* all omega */
        output_90 << data_phi_90[i][j] << " \t";
      output_90 << std::endl;
    }



  /* verbose information for the user: program is finished */
  std::cout << "Done" << std::endl;
  return 0;
}


