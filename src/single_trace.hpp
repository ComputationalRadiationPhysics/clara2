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
#include <string>
#include <cassert>
#include <cstdlib>



#ifndef SINGLE_TRACE_RPAUSCH
#define SINGLE_TRACE_RPAUSCH

using namespace std;

#include "vector.hpp"

#include "detector_e_field.hpp"
#include "detector_dft.hpp"
#include "detector_fft.hpp"

#include "import_from_file.hpp"

/* only needed for near field calculation 
 * REMOVE ? */
/* #include"large_index_storage.hpp" */

#include "physics_units.hpp"


/**
 * calculates a single spectra for only one trace and one direction
 *
 * @param one_line pointer to trajectory data
 * @param linenumber number of data points
 * @param all_omega pointer to frequency values
 * @param all_spectrum pointer to memory for spectra
 * @param N_all_spec maximum number of spectral data allocated
 * @param theta_offset offset of angle theta (used to set angle)
 * @param phi_offset offset of angle phi (used to set angle)
 **/
int single_trace(const one_line* data, 
                 const unsigned int linenumber,
                 const double* all_omega, 
                 double* all_spectrum, 
                 const unsigned N_all_spec,
                 const double theta_offset = 0.0, 
                 const double phi_offset = 0.0);



/* TO DO: this should be in a separate file - ISSUE #15 */
/**
 * check whether a file exists or not 
 *
 * @param filename pointer to array containing file location
 * @return Returs true if file exists, otherwise false.
 **/
bool file_exists(const char *filename);

#endif
