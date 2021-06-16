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


#pragma once

#include "vector.hpp"


/** \brief class for a point-like user-defined observation detector storing the signal externally 
 * Calculate E_adv at advanced time t_adv
*/
class Detector_uop
{
public:
  /** constructor for a point-like detector 
    * to calculate the electrical field on user-defined observation * plane
    *
    *  @param obs_vector   = vector in observation direction, 
    *  Vec(r_obs)
    *  @param N_data   = number of time step in trace
    *  @param t_obs_length = number of time step for user-defined *  observation plane
    */
  Detector_uop(const R_vec obs_vector, 
               const unsigned N_data
               );

  /** destructor */
  ~Detector_uop();

  /** method to calculate electrical field for one time step
    *
    * @param r = position
    * @param beta = speed / speed_of_light
    * @param dot_beta = time derivative of beta
    * @param t_part = time
    * @param delta_t = time step width
    */
  void calc_eField(const R_vec r_0,
                       const R_vec beta_0,
                       const R_vec dot_beta_0,
                       const double gamma,
                       const double t_part_0,
                       const double delta_t);
 
 /* Functions for interpolation */
 /** get N_data */
  unsigned get_number_data();

/** get time(int numberOfStep) */
  double get_timestep(unsigned a);

/* get data[a] */
  R_vec get_data(unsigned a);

//data:
private:

  const R_vec obs_vector; /* observation direction */
  const unsigned N_data; /* number of time steps in trace */

  const double t_obs_start; /* starting time for observation in detector plane */
  const double t_obs_end; /* ending time for observation in detector plane */
  unsigned counter; /* number of time steps analyzed */

  const double factor; /* SI unit factor for calculating electrical field */

  double* time; /* pointer to advanced time, t_adv */
  R_vec* data; /* pointer to electrical field at advanced times */

};
