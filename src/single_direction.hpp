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
