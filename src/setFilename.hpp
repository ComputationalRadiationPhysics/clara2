/**
 * Copyright 2016 Richard Pausch
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

#include "settings.hpp"


void setFilename(char* filename, 
                 const char* templateString,
                 const unsigned int trace_id, 
                 const unsigned int N_char_filename)
{
  if(sprintf(filename,
             templateString,
             trace_id)
     > int(N_char_filename))
  {
    /* throw warning when buffer is to small for path name */
    std::cerr << "filename buffer too small!!! " << std::endl;
    throw "Buffer to small!";
  }
}

