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


#include "detector_uop.hpp"

#include "../settings.hpp"
#include "physics_units.hpp"
#include "utilities.hpp"
#include <stdio.h>



// Constructor and Destructor:

/** constructor for a point-like detector using a Fast Fourier Transform
  * to calculate the radiation spectra
  *
  * @param obs_vector   = vector in observation direction, 
  *  Vec(r_obs)
  *  @param N_data   = number of time step in trace
  */
Detector_uop::Detector_uop(const R_vec obs_vector,
                           const unsigned N_data)
  : obs_vector(obs_vector),
    N_data(N_data-3),
    t_obs_start(param::t_start),
    t_obs_end(param::t_start+param::dt_obs*param::Nt_obs),
    counter(0),
    factor(-1.0*phy::q / (4.* M_PI *phy::epsilon_0)),
    time(0),
    data(0)
{
  time = new double[N_data]; /* advanced time */
  data = new R_vec[N_data]; /* electrical field at advanced times */
}

/** destructor */
Detector_uop::~Detector_uop()
{
  delete[] time;
  delete[] data;
}


// add to spectrum methods:

/** method to calculate electrical field for one time step
  *
  * @param r = electron position
  * @param beta = speed / speed_of_light
  * @param dot_beta = time derivative of beta
  * @param t_part = time
  * @param delta_t = time step width
  */
void Detector_uop::calc_eField(const R_vec r,
                               const R_vec beta,
                               const R_vec dot_beta,
                               const double gamma,
                               const double t_part,
                               const double delta_t)
{
  /* calculate electric field vector */
  const R_vec obs_direction = obs_vector-r;
  const R_vec n_unit = obs_direction.unit_vec();
  const double obs_direction_mag = obs_direction.mag();

  const R_vec fou1 = (n_unit-beta) / (util::square(obs_direction_mag*gamma) * util::cube(1. - beta * n_unit));
  const R_vec fou2 = (n_unit%((n_unit-beta)%dot_beta))
                     / (util::cube(1. - beta * n_unit) 
                     * phy::c *obs_direction_mag);

  /* advanced time */
  time[counter] = t_part + obs_direction_mag/phy::c;
  //printf("for timestep of electron: %.4e the observation vecotr amplitude is %3.8f \n", t_part, obs_direction_mag);


  /* abort if advanced time is not a subset of observation time */
  assert((time[counter] > t_obs_start) && "Starting advanced time of trace data is smaller than starting observation time");
  assert((time[counter] < t_obs_end) && "Ending advanced time of trace data is bigger than ending observation time"); 

  /* real vector amplitude */
  data[counter] = factor*(fou2); // + fou2);
  //printf("for timestep: %.4e: the eFields are: %3.8f, %3.8f,  %3.8f \n", time[counter], data[counter][0], data[counter][1], data[counter][2]);

  ++counter; /* count processed time steps */
}

unsigned Detector_uop::get_number_data()
{
  return N_data;
}

double Detector_uop::get_timestep(unsigned a)
{
  return time[a];
}

/* get data[a] */
R_vec Detector_uop::get_data(unsigned a)
{
  return data[a];
}
