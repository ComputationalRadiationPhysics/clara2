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
//#include <sstream>


#ifndef LOAD_TXT_RPAUSCH
#define LOAD_TXT_RPAUSCH

#include "import_from_file.hpp"


using namespace std;



//! \brief function returning the number of lines of a file
/*! @param target a 'const char*' containing file location */
unsigned linecounter(const char* target) 
{
    unsigned counter = 0; // line counter
    std::ifstream file(target); // target file
    while (!file.eof()) {
        if (file.get() == '\n') { ++counter; }
    }
    file.close(); // unnecessary but OK
    return counter + 1; // add one for last line without '\n'
}





void load_txt( const char target[], const unsigned linenumber, one_line* data)
{  
    ifstream file(target); // file to load
    
    string storage; // storage of short strings
    for (unsigned i=0; !file.eof(); ++i) 
      {
	// check for FORTRAN error: 1.234e-123 is stored wrongly as 1.234-123 !
#ifndef CHECK_FOR_FORTRAN_ERROR                
	file >> data[i/7].intern_data[i%7];
#else
	file >> storage;
        // going through the string and checking if the Mathematica output
        // error occurred:
        for(unsigned j=1; storage[j] != '\0'; ++j){ // ignoring a first sign
	  if (storage[j] == '-' && storage[j-1] != 'e'){
	    std::string str1, str2;
	    str1 = storage.substr(0, j); // string before error
	    str2 = storage.substr(j);    // string after error
	    storage = str1 + 'e' + str2; // corrected output
	  }
        }
        data[i/7].intern_data[i%7] = atof(storage.data());
        // storing data (7 doubles per line) (string to  double)
#endif
    }
    
    
    file.close();
    //std::cout << "Done\n\n";

    /////////////////////////////////////////////////////////////
    // GOT DATA /////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////
}


#endif

