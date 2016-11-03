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

//! \brief class for a point-like detector storing the signal externally
class Detector_dft
{
public:
    //! \brief constructor for a spectral detector
    /*!
       @param n_unit   = unit vector in direction of energy deposition
       @param delta_t  = time step of odint
     */
  Detector_dft(const R_vec n_unit,
               const unsigned spek_length,
               const double omega_max);

  Detector_dft(const R_vec n_unit,
               const unsigned spek_length,
               const double* omega);


  ~Detector_dft();

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

  double get_spectrum(unsigned a,
                      unsigned b);

  double energy();

// debugging
  Vector<std::complex<double>, 3>* spektrum;

private:
  // data:
  const R_vec n_unit;
  const  unsigned spek_length;

  double* spektrum_mag;
  double* frequency;

  // set frequency methods:
  inline void set_frequency(const double omega_max);
  inline void set_frequency(const double* omega);

};
