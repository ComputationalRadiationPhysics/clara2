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
 */




#include <iostream>
#include <string>
#include <cstdlib>
#include <fstream>
#include <cmath>


/* TO DO: this should be in a separate file - ISSUE #15 */
/**
 * This function checks if a file exists on the hard drive.
 *
 * @param filename string containg the path and filename to bechecked
 * @return Returns true if file exists, otherwise false.
 */
bool file_exists(const char *filename)
{
  std::ifstream infile(filename);
  return infile;
}






int main()
{

  std::cout << "Convert single file data to matrix: " << std::flush;

  // ----------- set parameters ------------ //
  
  const unsigned N_omega = 2048; // could be read from linenumber input file
  double omega[N_omega];


  const unsigned N_theta = 120; // number of different theta angles
  double theta[N_theta];
  double theta_max =  1.14594939 /** 120*/ ; // degree = 1.0/gamma
  for(unsigned i=0; i<N_theta; ++i)
    {
      theta[i] = (double)i / N_theta * theta_max;
    }

  const unsigned N_phi = 2; // number of different phi angles
  double phi[N_phi] = {0.0, 90.0};




  double data_phi_0[N_theta][N_omega]; // matrix data  for phi=0 degrees
  double data_phi_90[N_theta][N_omega];// matrix data  for phi=90 degrees



  // -------- end set parameters -------- //




  // -------- begin load data -------- //

  char buffer[10];
  for(unsigned i=0; i< N_phi*N_theta; ++i)
    {
      sprintf(buffer, "%3.5f", theta[i%N_theta]);
      std::string theta_st = buffer; // string with correct format for theta
      sprintf(buffer, "%3.5f", phi[i/N_theta]);
      std::string phi_st = buffer; // string with correct format for phi
      std::string filename = "my_output_theta=" + theta_st + "_phi=" + phi_st + ".dat";
      if(!(file_exists(filename.c_str()))) // check if file exists
	std::cout << "  ERROR: file " << filename << " does not exists" << std::endl;
      else // file exists:
	{
	  std::ifstream data(filename.c_str()); // load file
	  std::cout << "loading file: " << filename << std::endl;
	  if(i==0) // for fist file save also omega values
	    {
	      for(unsigned j=0; j<N_omega; ++j)
		data >> omega[j] >> data_phi_0[i][j];
	    }
	  else if(i<N_theta) // for rest of phi=0
	    {
	      double dummy;
	      for(unsigned j=0; j<N_omega; ++j)
		{
		  data >> dummy >> data_phi_0[i][j];
		  if(dummy != omega[j]) // check if frequencies are equal
		    std::cout << "  ERROR: frequency discrepancy: " << filename
			      << " - " << omega[j] << " != " << dummy << std::endl;
		}
	    }
	  else // all files for phi=90 --> same procedure
	    {
	      double dummy;
	      for(unsigned j=0; j<N_omega; ++j)
		{
		  data >> dummy >> data_phi_90[i%N_theta][j];
		  if(dummy != omega[j])
		    std::cout << "  ERROR: frequency discrepancy: " << filename
			      << " - " << omega[j] << " != " << dummy << std::endl;
		}
	    }
	  data.close();
	}
    } 


  // -------- end load data --------- //



  std::ofstream output_0("matrix_phi_0.dat");
  for(unsigned i=0; i<N_theta; ++i)
    {
      for (unsigned j=0; j<N_omega; ++j)
	output_0 << data_phi_0[i][j] << " \t";
      output_0 << std::endl;
    }


  std::ofstream output_90("matrix_phi_90.dat");
  for(unsigned i=0; i<N_theta; ++i)
    {
      for (unsigned j=0; j<N_omega; ++j)
	output_90 << data_phi_90[i][j] << " \t";
      output_90 << std::endl;
    }








  
  std::cout << "Done" << std::endl;
  return 0;
}


