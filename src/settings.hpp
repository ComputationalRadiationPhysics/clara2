/**
 * Copyright 2016 Alexander Koehler, Richard Pausch
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

namespace param
{
const double omega_max                 = 3.0e19;      /* maximum of plotted frequency Hz */
const double theta_max                 = 1.14594939;  /* maximum of theta in degree */
const unsigned int N_spectrum          = 2048; /* number of frequencies "omega" */
const unsigned int N_theta             = 120;     /* number of directions in first angle "theta" */
const unsigned int N_phi               = 2;         /* number of directions in second angle "phi" */
const unsigned int N_trace             = 1;    /* maximum number of traces */

const unsigned int fft_length_factor   = 1; /* needs to be a power of two */

const bool ascii_output = true; /* output spectra as ascii text */

const unsigned int N_char_filename=256; // number of characters
/* template for input traces (using C formatting): */
const char traceFileTemplate[] = "trace_%04d.txt";
/* template for output spectra for each trace (using C formatting): */
const char outputFileTemplate[] = "my_spectrum_trace%06d.dat";
const char outputFileTemplate_uop[] = "my_eField_trace%06d.dat";

///////////////////////////////////////////////////////////////
///////////////////////// by PENGHAO///////////////////////////

/* Parameters for User-defined Observation Plane electrical field calculation */

const bool USE_uop                     = true; /* Whether to use User-defined Observation Plane electrical field calculation*/ 

const double X_uop                     = 0.0;    /* X position of user defined observation plane, Unit:[m]. */

/* Box length of observation plane, Unit:[m]. Note that coordinates of observation point are associated with observation angle by: (y, z) = R*(sin(theta)*sin(phi), sin(theta)*cos(phi)), where R=|X_part-X_uop|, X_part is x coordinate of electron.*/
const double Y_max                     = 30.00625;   // Last point is not reached, the actual Y_max is 2e-6     

const double Y_min                     = -0;

const double Z_max                     = 0.0;     // // Last point is not reached, the actual Z_max is 2e-6

const double Z_min                     = 0.0;

/* Resolution of the observation plane */
const unsigned int Ny                  = 4801;

const unsigned int Nz                  = 1;

const double t_start                   = -1.5e-7; /* When to start recording electrical field on observation plane, usually the time first photon arrives at observation plane, Unit:[s].*/
const double dt_obs                    = 1e-10; /* Timestep to record electrical field, Unit:[s]. This is associated with PI/omega_max */

const unsigned Nt_obs                  = 3000;   /*number of record timesteps. Notice users need to make sure t_start<t_part_start, t_end>t_part_end, where t_end= t_start+Nt_obs*dt_obs, t_part_start and t_part_end are arriving time of first photon and last photon emitted by a single electron at observation plane .  

For case of inverse thomson scattering, emission time of a single electron is T_emi=N_osc*lambda_x/c, N_osc is number of electron oscillations, lambda_x is wavelength of emitted x-ray, which is ususally on order of attosecond. Emission time of an electron beam is mainly determined by bunch time length T_bunch, so users need to make sure Nt_obs*dt_obs>T_bunch+2*T_emi.
*/
const double t_pick =  0.0;

const unsigned Nt_pick = 1500;
///////////////////////// by PENGHAO///////////////////////////
///////////////////////////////////////////////////////////////
namespace input
{
  const unsigned int index_time = 6; /* column id of time */

  const unsigned int index_pos_x = 0; /* column id of x */
  const unsigned int index_pos_y = 1; /* column id of y */
  const unsigned int index_pos_z = 2; /* column id of z */

  const unsigned int index_beta_x = 3; /* column id of beta_x */
  const unsigned int index_beta_y = 4; /* column id of beta_y */
  const unsigned int index_beta_z = 5; /* column id of beta_z */

  /* unit conversion to SI units */
  const double convert_time = 1.0; /* s -> s */
  const double convert_pos = 1.0; /* m -> m */
  const double convert_beta = 1.0; /* None -> None */
}

// explicit for process_data
const unsigned int N_omega             = N_spectrum; /* number of frequencies */
const unsigned int index_files_first   = 0; /* start index for reading trajectory spectra */
const unsigned int index_files_last    = N_trace; /* final index for reading trajectory spectra */
/* template for output total (incoherent) spectra  (using C formatting for phi index): */
const char output_pattern[] = "my_spectrum_all_%03d.dat";
const char output_pattern_uop[] = "my_efield_all_%03d.dat";

} // end namespace param
