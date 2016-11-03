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

#include "fileExists.hpp"

/**
 * check whether a file exists or not
 *
 * @param filename pointer to array containing file location
 * @return Returs true if file exists, otherwise false.
 **/
bool file_exists(const char *filename)
{
  std::ifstream ifile(filename);
  return ifile;
}
