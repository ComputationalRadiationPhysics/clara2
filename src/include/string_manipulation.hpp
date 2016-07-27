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
#include <cstdlib>
#include <cstdio>

#ifndef formating_rpausch
#define formating_rpausch

void formating(std::string& index_s, const std::string& digits_s)
{
 // ugly formatting:
  unsigned index_u = atoi(index_s.c_str());  // make a number
  unsigned digits_u = atoi(digits_s.c_str());// make a number
  
  unsigned dummy=1;
  for(unsigned i = 0; i< digits_u; ++i)
    dummy *=10;

  if(index_u >= dummy) std::cout << " WARNING: index too large" << std::endl;

  std::string format = ("%0" + digits_s + "d"); // format style
  char* index_c = new char[digits_u+5]; // create char dummy
  sprintf(index_c, format.c_str(), index_u); // convert from unsigned to char
  index_s = index_c; // convert char to string
  delete[] index_c; // delete char dummy
  // end formatting
}


double get_double(std::string txt)
{
  return atof(txt.c_str());
}



unsigned get_unsigned(std::string txt)
{
  return (unsigned)atoi(txt.c_str());
}

#endif
