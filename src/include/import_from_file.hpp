/**
 * Copyright 2014 Richard Pausch
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




#include <fstream>

#ifndef IMPORT_FROM_FILE_RPAUSCH
#define IMPORT_FROM_FILE_RPAUSCH


//! \brief simple container to store data from the Clara trace
/*! usage: one_line x[number of data lines]; then x[i].intern_data[0-6] */
struct one_line {
    double intern_data[7]; // simple data structur two handel 7 doubles per line
};

#endif

