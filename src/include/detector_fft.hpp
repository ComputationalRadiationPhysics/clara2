/**
 * Copyright 2014-2018 Richard Pausch
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


#pragma once

#include "vector.hpp"


/** \brief class for a point-like detector storing the signal externally */
class Detector_fft
{
public:
  /** constructor for a point-like detector using a Fast Fourier Transform
    * to calculate the radiation spectra
    *
    *  @param n_unit   = unit vector in direction of energy deposition
    *  @param N_data   = number of time step in trace
    */
  Detector_fft(const R_vec n_unit, 
               const unsigned N_data);

  /** destructor */
  ~Detector_fft();

  /** method to add an radiation amplitude for one time step
    *
    * @param r = position
    * @param beta = speed / speed_of_light
    * @param dot_beta = time derivative of beta
    * @param t_part = time
    * @param delta_t = time step width
    */
  void add_to_spectrum(const R_vec r_0,
                       const R_vec beta_0,
                       const R_vec dot_beta_0,
                       const double t_part_0,
                       const double delta_t);


  /* ISSUE # 74 - interface too crowded */
  /** method to add an radiation amplitude for one time step
    *
    * @param r = position
    * @param p = momentum
    * @param dot_p = time derivative of p
    * @param beta = speed / speed_of_light
    * @param gamma = relativistic gamma factor
    * @param dot_gamma = time derivative of gamma
    * @param t_part = time
    * @param delta_t = time step width
    */
  void add_to_spectrum(const R_vec r_0,
                       const R_vec p_0,
                       const R_vec dot_p_0,
                       const R_vec beta_0,
                       const double gamma_0,
                       const double dot_gamma_0,
                       const double t_part_0,
                       const double delta_t);

  /** calculate the spectrum from all added amplitudes */
  void calc_spectrum();

  /* ISSUE #75 - separate get spectrum and get frequency */
  /** get spectrum or frequency
    *
    * @param a = id between [0, number of frequencies]
    * @param b = int return 0->frequency, 1->spectra in Js
    */
  double get_spectrum(unsigned a, unsigned b);

  /** return total energy (in calculated spectral range) */
  double energy();

  /** return the number of frequencies till the Nyquist frequency
    * (number of half of all frequency bins)
    */
  unsigned half_frequency();

//data:
private:

  const R_vec n_unit; /* observation direction */
  const unsigned N_data; /* number of time steps in trace */

  unsigned spec_length; /* number of frequencies for spectra */
  unsigned counter; /* number of time steps analyzed */

  double* time; /* pointer to retarded time */
  R_vec* data; /* pointer to real vector amplitudes */

  double* spectrum_mag; /* pointer to absolute spectra */

public:
  double* frequency; /* pointer to frequencies */

};
