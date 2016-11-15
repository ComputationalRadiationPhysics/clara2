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



#include <iostream>
#include <cassert>
#include <complex>
#include <cmath>

#include "vector.hpp"
#include "physics_units.hpp"
#include "utilities.hpp"

#pragma once

/** Detector class that computes the spectra via a
  * (non-equidistat) discrete Fourier transform
  */
class Detector_dft
{
public:
  /** constructor for a spectral detector
    *
    * @param n_unit      = unit vector pointing in observation direction
    * @param spec_length = number of frequenies to compute
    * @param omega_max   = maximum frequency -> frequency range [0, omega_max]
   */
  Detector_dft(const R_vec n_unit,
               const unsigned spec_length,
               const double omega_max);

  /** constructor for a spectral detector
    *
    * @param n_unit      = unit vector pointing in observation direction
    * @param spec_length = number of frequenies to compute
    * @param omega       = pointer to frequencies to calculate
   */
  Detector_dft(const R_vec n_unit,
               const unsigned spec_length,
               const double* omega);


  /** destructor */
  ~Detector_dft();

  /** method to add an radiation amplitude for one time step
    *
    * @param r_0 = position
    * @param beta_0 = speed / speed_of_light
    * @param dot_beta_0 = time derivative of beta_0
    * @param t_part_0 = time
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
    * @param r_0 = position
    * @param p_0 = momentum
    * @param dot_p_0 = time derivative of p_0
    * @param beta_0 = speed / speed_of_light
    * @param gamma_0 = relativistic gamma factor
    * @param dot_gamma_0 = time derivative of gamma_0
    * @param t_part_0 = time
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


  /** get spectrum or frequency
    *
    * @param a = id between [0, number of frequencies]
    * @param b = int return 0->frequency, 1->spectra in Js
    */
  double get_spectrum(unsigned a,
                      unsigned b);

  /** return total energy (in calculated spectral range) */
  double energy();


private:
  // data:
  const R_vec n_unit; /* observation direction */
  const unsigned spec_length; /* number of frequencies to calculate */
  Vector<std::complex<double>, 3>* spectrum; /* complex amplitudes */
  double* spectrum_mag; /* resulting spectrum */
  double* frequency; /* pointer to frequency values */

  // set frequency methods:

  /** set frequency based on maximum frequency
   * to cover frequency range [0, omega_max]
   *
   * @ param omega_max = maximum frequency
   */
  inline void set_frequency(const double omega_max);

  /** set frequency using pointer to frequency values
   *
   * @ param omega = double pointer to frequency values
   */
  inline void set_frequency(const double* omega);

};
