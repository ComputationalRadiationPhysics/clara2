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

/**
 * function that calculates spectra in different directions for
 * a single particle trace
 *
 * @param trace_id a unique id which which the trajectory file
 *                  can be identified
 * @return error code
 **/
int all_directions(const unsigned int trace_id);

/////////////////// modified by PENGHAO////////////////////
/**
 * function that calculates electrical field on the user-defined
 * observation plane in different directions for a single particle
 *  trace
 *
 * @param trace_id a unique id which identifies the trajectory file
 *
 * @return error code
 **/
int all_directions_uop(const unsigned int trace_id);
