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


#ifndef PARALLELJOBSRPAUSCH
#define PARALLELJOBSRPAUSCH

#if __PARALLEL_SETTING__ == 1
#include "mpi.h"
#endif


int start_array(int* numtasks, int* rank)
{
  if(!numtasks)
    return 2;

  if(!rank)
    return 3;

#if __PARALLEL_SETTING__ == 1
  // MPI used for parallization
  int rc = MPI_Init(NULL, NULL);

  if (rc != MPI_SUCCESS)
    {
      printf("Error starting MPI program. Terminating program!\n");
      MPI_Abort(MPI_COMM_WORLD, rc);
      return 1;
    }

  MPI_Comm_size(MPI_COMM_WORLD, numtasks);
  MPI_Comm_rank(MPI_COMM_WORLD, rank);

#elif __PARALLEL_SETTING__ == 2
  // Array jobs used for parallization
  char* dump;
  dump = getenv("PBS_ARRAYID");
  *rank = atoi(dump);
  dump = getenv("ARRAYMAX");
  *numtasks = atoi(dump);

#else
  // non of the above selected
  #error parallel setting not suported
#endif

  printf("Number of tasks= %d My rank= %d\n", *numtasks, *rank);

  return 0;
}

int end_array(void)
{
#if __PARALLEL_SETTING__ == 1
  // MPI used for parallization  
   MPI_Finalize();

#elif __PARALLEL_SETTING__ == 2
  // Array jobs used for parallization

#else
  // non of the above selected
  #error parallel setting not suported
#endif

   return 0;
}

/*
int check_break(void)
{
  char stop[] = "break.now";
  FILE* file = fopen(stop, "r");
  if(file != 0)
    {
      fclose(file);
      //perror("breakpoint found \n");
      return 1;
    }

  return 0;
}
*/

#endif

