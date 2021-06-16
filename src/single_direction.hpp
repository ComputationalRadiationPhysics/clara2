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

#include "import_from_file.hpp"
#include "vector.hpp"


/**
 * calculates a single spectra for only one trace and one direction
 *
 * @param data pointer to trajectory data
 * @param linenumber number of data points
 * @param all_omega pointer to frequency values
 * @param all_spectrum pointer to memory for spectra
 * @param N_all_spec maximum number of spectral data allocated
 * @param theta_offset offset of angle theta (used to set angle)
 * @param phi_offset offset of angle phi (used to set angle)
 **/

int single_direction(const one_line* data,
                     const unsigned int linenumber,
                     const double* all_omega,
                     double* all_spectrum,
                     const unsigned N_all_spec,
                     const double theta_offset = 0.0,
                     const double phi_offset = 0.0);


/**
 * calculates a single electrical filed for only one trace and one direction
 *
 * @param data pointer to trajectory data
 * @param linenumber number of data points
 * @param all_t_obs pointer to time values
 * @param all_eField pointer to memory for eField
 * @param N_all_eField maximum number of observation time allocated
 * @param x_offset x offset of observation plane (used to set
 * observation direction)
 * @param y_offset y offset of observation plane (used to set
 * observation direction)
 * @param z_offset z offset of observation plane (used to set
 * observation direction)
 **/
int single_direction_uop(const one_line* data,
                     const unsigned int linenumber,
                     const double* all_t_obs,
                     R_vec* all_eField,
                     const unsigned N_all_eField,
                     const double x_offset = 0.0, 
                     const double y_offset = 0.0,
                     const double z_offset = 0.0);
