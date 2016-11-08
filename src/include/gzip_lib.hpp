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



#include <zlib.h>
#include <stdio.h>
#include <iostream>

#pragma once


/* TO DO: not dealing with gzib files - ISSUE #19 */
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



/* TO DO: not dealing with gzib files - ISSUE #19 */
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



/* TO DO: is compression really used in the code? - ISSUE #5 */
/**
 * write data to file using compression
 *
 * @param data pointer to data that should be stored on disk
 * @param size number of bytes that are in data and should be
 *             put in the file
 * @param filename pointer to char-array containing path and
 *                 filename of the file to put data in
 * @param number "unsigned int" that  is written to the
 *               beginning of the file
 * @return returns zero in case of success, no error handling
 **/
int compress_data(void* data,
                  long unsigned int size,
                  char* filename,
                  unsigned int number)
{
  /* create compressed file handler (write binary="wb") */
  gzFile outfile = gzopen(filename, "wb");

  /* write number in front of file */
  gzwrite(outfile, &number, sizeof(unsigned int));
  /* write data to file */
  gzwrite(outfile, data, size);

  /* close file handler */
  gzclose(outfile);

  return 0;
}



/* TO DO: is compression really used in the code? - ISSUE #5 */
/**
 * append data to a compressed file
 *
 * @param data pointer to data that should be stored on disk
 * @param size number of bytes that are in data and should be
 *             put in the file
 * @param filename pointer to char-array containing path and
 *                 filename of the file to append data to
 * @return returns zero in case of success, no error handling
 **/
int compress_data_append(void* data,
                         long unsigned int size,
                         char* filename)
{
  /* create compressed file handler (append binary="ab") */
  gzFile outfile = gzopen(filename, "ab");

  /* append data to file */
  gzwrite(outfile, data, size);

  /* close file handler */
  gzclose(outfile);

  return 0;
}



/* TO DO: is this function still necessary? - ISSUE #5 */
/**
 * specific function to create a gzFile to access
 * data from the "bubbleStressTest" simulation
 *
 * @param pFile empty "gzFile" that should point
 *              to simulation data
 * @param od_or_even select whether the file lies in
 *                   bigOutput1 (odd) or bigOutput2
 *                   (even)
 * @param time_id simulation time step to identify
 *                file
 * @param core_id CPU-core number for file name
 **/
void create_gzFile(gzFile& pFile,
                   int od_or_even,
                   unsigned int time_id,
                   unsigned int  core_id)
{
  /* set up data file: */

  /* create file name from "bubbleStressTest" path*/
  char filename[256];
  if(sprintf(filename, /* char-array to put path into */
             "/net/cns/projects/bubbleStressTest/bigOutput%1d/e%05d_%03d.dat.gz",
             od_or_even, /* select path depending on odd or even */
             time_id, /* identify file by time step */
             core_id /* CPU core identification used for file naming */
             ) > 254)
   {
     /* throw error in case the buffer is to small */
     std::cerr << "buffer  to small!!! " << std::endl;
     /* TO DO: error handling - ISSUE #20 */
     throw "Buffer to small!";
   }

  /* create gzFile handler for read-only */
  pFile = gzopen(filename, "rb");

  /* verify if file was created correctly */
  if(pFile == NULL) /* if error occurred */
  {
    pFile = gzopen(filename, "rb"); /* try again */
    if(pFile == NULL) // try again
    {
      /* if second attempt fails too throw error */
      std::cerr << "Could not open file" << std::endl;
      /* TO DO: error handling - ISSUE #20 */
      throw "no file found";
    }
  }
}
