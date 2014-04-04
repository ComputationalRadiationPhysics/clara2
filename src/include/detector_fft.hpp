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




#include <iostream>
#include <cassert>
#include <complex>
#include <cmath>

#include "vector.hpp"
#include "physics_units.hpp"
#include "utilities.hpp"
#include "fft_ned.hpp"



#ifndef DETECTOR_FFT_RPAUSCH
#define DETECTOR_FFT_RPAUSCH

//! \brief class for a point-like detector storing the signal externaly
class Detector_fft
{
public:
  //! \brief constructor for a point-like detector
  /*! 
     @param n_unit   = unit vector in direction of energy deposition
     @param delta_t  = timestep of odint
  */
  Detector_fft(R_vec n_unit, unsigned N_data);

  ~Detector_fft();

  void add_to_spectrum(const R_vec r_0, 
		       const R_vec beta_0,
		       const R_vec dot_beta_0,
		       const double t_part_0,
		       const double delta_t);


  void add_to_spectrum(const R_vec r_0, 
		       const R_vec p_0, 
		       const R_vec dot_p_0,
		       const R_vec beta_0,
		       const double gamma_0,
		       const double dot_gamma_0,
		       const double t_part_0,
		       const double delta_t);
	

  void calc_spectrum();

	
  double get_spectrum(unsigned a, unsigned b);

  double energy();

  unsigned half_frequency();

private:

  //data:
  
  const R_vec n_unit;
  // delta_t; // neccesairy for integration
  const unsigned N_data;

  //public:  // --> better to private !!!
  unsigned spek_length;
  unsigned counter;

  double* time;
  R_vec* data;

  double* spektrum_mag;

public:
  double* frequency;
	
};


#endif

