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



#include <zlib.h>
#include <stdio.h>
#include <iostream>

#ifndef GZIP_LIB_RPAUSCH
#define GZIP_LIB_RPAUSCH




int store_data(void* data, long unsigned int size, char* filename)
{
  FILE* outfile = fopen(filename, "wb");

  if(outfile == NULL)
    {
      std::cerr << "could not create file: " << filename << std::endl;
      throw "could_not_create_file";
    }

  fwrite(data, 1, size, outfile);
  fclose(outfile);
  
  return 0;
}


int read_data(void* data, long unsigned int size, char* filename)
{
  FILE* outfile = fopen(filename, "rb");

  if(outfile == NULL)
    {
      return 1;
    }

  fread(data, 1, size, outfile);
  fclose(outfile);
  
  return 0;
}






int compress_data(void* data, long unsigned int size, char* filename, 
    unsigned int number)
{
  gzFile outfile = gzopen(filename, "wb");
  gzwrite(outfile, &number, sizeof(unsigned int));
  gzwrite(outfile, data, size);
  gzclose(outfile);
  return 0;
}



int compress_data_append(void* data, long unsigned int size, char* filename)
{
  gzFile outfile = gzopen(filename, "ab");
  gzwrite(outfile, data, size);
  gzclose(outfile);
  return 0;
}


void create_gzFile(gzFile& pFile, int od_or_even, unsigned int time_id, 
  unsigned int  core_id)
{
  /// set up data file:
  char filename[256];
  if(sprintf(filename, "/net/cns/projects/bubbleStressTest/bigOutput%1d/e%05d_%03d.dat.gz", 
    od_or_even, time_id, core_id) > 254)
    {
      std::cerr << "buffer  to small!!! " << std::endl;
      throw "Buffer to small!";
    }

  //printf( "load file: \t%s \n", filename);

  pFile = gzopen(filename, "rb"); // location, Read-only

  if(pFile == NULL)
    {
      pFile = gzopen(filename, "rb");
        if(pFile == NULL) // try again
          {
            std::cerr << "Could not open file" << std::endl;
            throw "no file found";
          }
    }
}





#endif
