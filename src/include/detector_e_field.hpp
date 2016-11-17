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
#include "large_index_storage.hpp"
#include "physics_units.hpp"
#include "utilities.hpp"


#pragma once


/** \brief class for a point-like detector storing the electric field (signal) externally */
class Detector_e_field
{
public:
  /** constructor for electric field detector (at position in space)
    *
    * @param detector  = R_vec with observation position
    * @param delta_t   = time step width
    * @param N_sig     = maximum number of field entries to the detector
    * @param start_sig = start time at which the electric field should
    *                    be "recorded"
    */
  inline Detector_e_field(R_vec detector,
                          double delta_t,
                          unsigned N_sig,
                          unsigned start_sig );

  /** compute integer index in which the retarded electric field (signal)
    * should be written to
    *
    * @param t_signal = time at which the electric field arrives at
    *                   detector position
    * @return unsigned int index in which the time fits best in
    *         the signal class
    */
  inline unsigned delay_index(double t_signal);

  /** compute the time at which the electric field emitted by
    * a charged particle at position r and time t arrives at the
    * detector position (retarded time)
    *
    * @param r = R_vec position of particle
    * @param t = double current time (of the particle)
    * @return signal arrival time (retarded time)
    */
  inline double t_signal(R_vec r, double t);

  /* ISSUE #74 - clean up interface */
  /** compute electric field value(s) at detector position
    * and store the field values using two time steps
    * (subscript 0 and 1)
    *
    * @param r_0 = particle position at t_part_0
    * @param r_1 = particle position at t_part_0 + delta_t
    * @param p_0 = particle momentum at t_part_0
    * @param p_1 = particle momentum at t_part_0 + delta_t
    * @param dot_p_0 = particle change in momentum at t_part_0
    * @param dot_p_1 = particle change in momentum at t_part_0 + delta_t
    * @param beta_0 = particle beta (v/c) at t_part_0
    * @param beta_1 = particle beta (v/c) at t_part_0 + delta_t
    * @param gamma_0 = particle relativistic gamma at t_part_0
    * @param gamma_1 = particle relativistic gamma at t_part_0 + delta_t
    * @param dot_gamma_0 = particle change in gamma at t_part_0
    * @param dot_gamma_1 = particle change in gamma at t_part_0 + delta_t
    * @param t_part_0 = time at 'first' time step
    */
  void place(const R_vec r_0,
             const R_vec r_1,
             const R_vec p_0,
             const R_vec p_1,
             const R_vec dot_p_0,
             const R_vec dot_p_1,
             const R_vec beta_0,
             const R_vec beta_1,
             const double gamma_0,
             const double gamma_1,
             const double dot_gamma_0,
             const double dot_gamma_1,
             const double t_part_0);

  /** returns the counter of double counts minus no counts
    *
    * @return value of index counter
    */
  inline int count();


  /* data: */
  const double delta_t;  /* length of time steps */
  Large_index_storage<R_vec> signal;  /* container for electric field at detector */

private:
  /* data: */
  R_vec detector; /* location of the detector */
  int counter; /* the counter of double counts minus no counts
                * (info for user to optimze detector) not needed for calulations */

  /* methods: */

  /** simple interpolation between two values (f.e. location, speed)
    * at two time step with time difference delta_t
    *
    * @param x_0   = value at start point
    * @param x_1   = value at endpoint (one time step later)
    * @param t     = time as interpolation parameter (0 < t < delta_t)
    */
  template<typename V>
  V interpol(V r_0,
             V r_1,
             double t);

  /** calculates the electric field $\vec E$ based on Lienard Wiechert potential
    * @param e_R      = unit vector pointing from the electron to the detector
    * @param beta     = beta of the electron
    *                   $ \vec \beta = \frac{\vec v}{c} $
    * @param beta_dot = $ \operatorname{\frac{d}{dt}} \vec \beta $
    * @param gamma    = gamma value of the electron
    *                   $ gamma = \sqrt{\frac{1}{1-\vec \beta^2} } $
    * @param R       = distance between electron and detector
    * @return electric field at detector position
    */
  R_vec Lienard_Wiechert(const R_vec& e_R,
                         const R_vec& p,
                         const R_vec beta_dot_times_gamma,
                         const double& gamma,
                         const double& R);

};
