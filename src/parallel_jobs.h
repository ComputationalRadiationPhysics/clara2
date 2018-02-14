/**
 * Copyright 2014-2018 Richard Pausch
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


/* pre-compiler variable determines whether parallelization is done via
 * PBS (parallel batch system) on the cluster and used the cluster
 * scheduler or via MPI (message passing interface).
 * PBS is best use if you want to fill the cluster as best as possible
 * but as a overhead due to the time PBS needs to schedule each task.
 * If enough resources are available, MPI is faster and puts less load
 * an the PBS scheduler
 *
 * __PARALLEL_SETTING__ == 1 --> MPI
 * __PARALLEL_SETTING__ == 2 --> PBS array job
 */
#if __PARALLEL_SETTING__ == 1
/* only include mpi header only if MPI was selected */
#include "mpi.h"
#endif


/** this function emulates a PBS Array job task and either runs it
  * in a PBS array environment or in an MPI environment
  *
  * @param numtasks pointer to int where number of total parallel
  *                 tasks should be stored
  * @param rank pointer to int where the rank/task-id of a specific
  *             task should be stored
  * @return int with error code: 0 - >successful
  *                              1 -> MPI error
  *                              2 -> numtask not allocated
  *                              3 -> rank not allocated
  */
int start_array(int* numtasks,
                int* rank)
{
  /* check if numtask is a valid pointer */
  if(!numtasks)
    return 2;

  /* check if rank is a valid pointer */
  if(!rank)
    return 3;

/* in case MPI is used for parallelization */
#if __PARALLEL_SETTING__ == 1
  /* init MPI on this core */
  int rc = MPI_Init(NULL, NULL);

  /* check if MPI init was successful */
  if (rc != MPI_SUCCESS)
  {
    printf("Error starting MPI program. Terminating program!\n");
    MPI_Abort(MPI_COMM_WORLD, rc);
    return 1;
  }

  /* get MPI numtasks and rank */
  MPI_Comm_size(MPI_COMM_WORLD, numtasks);
  MPI_Comm_rank(MPI_COMM_WORLD, rank);

/*in case  PBS array jobs are used for parallelization */
#elif __PARALLEL_SETTING__ == 2
  char* dump; /* temporary memory pointer for reading bash
                 environment variables */
  /* get rank = PBS_ARRAYID */
  dump = getenv("PBS_ARRAYID");
  *rank = atoi(dump);
  /* get numtask = ARRAYMAX */
  dump = getenv("ARRAYMAX");
  *numtasks = atoi(dump);

/* throw compile time error if neither MPI nor PBS array jobs are selected */
#else
  #error parallel setting not suported
#endif

  printf("Number of tasks= %d My rank= %d\n", *numtasks, *rank);

  return 0;
}


/** function to clean up after parallel job (if needed)
  * @return int error code: 0 -> successful
  */
int end_array(void)
{
/* in case MPI is used for parallelization */
#if __PARALLEL_SETTING__ == 1
  MPI_Finalize();
/* in case PBS array jobs are used for paralellization */
#elif __PARALLEL_SETTING__ == 2
  /* nothing needs to be done */
/* throw compile time error if neither MPI nor PBS array jobs are selected */
#else
  #error parallel setting not suported
#endif

  return 0;
}


/** function to check if a specific file is found in the working
  * directory that forces a soft stop
  *
  * @return int 0 -> file not found (continue)
  *             1 -> file found (stop running)
  */
int check_break(void)
{
  char stop[] = "break.now"; /* file name that causes a stop */
  FILE* file = fopen(stop, "r"); /* try read access to file */
  if(file != 0) /* file found */
  {
    fclose(file); /* close file handler again */
    return 1; /* return 1 = file found */
  }

  return 0; /* file not found - return 0 */
}
