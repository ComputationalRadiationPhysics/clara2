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




#include <stdio.h>
#include <stdlib.h>

#include "parallel_jobs.h"
#include "all_directions.hpp"

int main(void)
{
  printf("start\n");

  int numtasks, rank;

  if(start_array(&numtasks, &rank) != 0)
    return 1;

  //////////////////////////////////////////////////
  //    do some work here
  //////////////////////////////////////////////////
  
  char * pHost;
  pHost = getenv("MYHOSTNAME");

  const char output_file[] = "my_output.txt-%d";
  const char error_file[] = "my_error.txt-%d";


  int return_value;
  char dump[1024] ={0};

  const unsigned int N_max = 2001;
  unsigned int i;
  for(i=rank; i<N_max; i += numtasks)
    {
      printf("this is job %5d of %5d jobs in the array (on %s = rank: %d)\n", i, N_max, pHost, rank);
      
      sprintf(dump, output_file, i);
      freopen(dump, "w", stdout);
      sprintf(dump, error_file, i);
      freopen(dump, "w", stderr);

      return_value = all_directions(i, "binary");
      
      fclose(stdout);
      fclose(stderr);
      freopen("/dev/tty", "w", stdout);
      freopen("/dev/tty", "w", stderr);



      /*
      if(check_break())
	{
	  printf("break point was set: i= %d  rank=%d  numtasks=%d \n", i, rank,  numtasks);
	  break;
	}
      */

    }



  ///////////////////////////////////////////////////
  //    work is done
  //////////////////////////////////////////////////
  
  end_array();
 

  return 0;
}
