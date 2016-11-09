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



#include <stdio.h>
#include <iostream>

#pragma once


/**
 * write any data as binary to file
 *
 * @param data pointer to data (array)
 * @param size size of the data in bytes
 * @param filename pointer to a char-array containing the filename
 *                 to store data in
 * @return gives zero on success, throws error if file can not
 *         be written to disk
 **/
int store_data(void* data,
               long unsigned int size,
               char* filename)
{
  /*  create file handler to write in binary format = "wb" */
  FILE* outfile = fopen(filename, "wb");

  if(outfile == NULL) /* in case file can not be written */
  {
    std::cerr << "could not create file: " << filename << std::endl;
    /* TO DO: error handling - ISSUE #20 */
    throw "could_not_create_file";
  }

  fwrite(data, 1, size, outfile);  /* write data to file */
  fclose(outfile); /* close file handler */

  return 0;
}



/**
 * read binary data from file
 *
 * @param data pointer to memory where data should be stored
 * @param size number of bytes to be read from file and put into data
 * @param filename pointer to char array with path and filename of
 *                 the file containing the data
 * @return zero on success, one in case the file can not be read
 **/
int read_data(void* data,
              long unsigned int size,
              char* filename)
{
  /*  create file handler to read in binary format = "rb" */
  FILE* outfile = fopen(filename, "rb");

  /* in case the file handler is erroneous (file does not exist,
   * etc.) return 1 */
  if(outfile == NULL)
  {
    /* TO DO: error handling - ISSUE #20 */
    return 1;
  }

  fread(data, 1, size, outfile);   /* read data from file*/
  fclose(outfile); /* close file handler */

  return 0;
}
