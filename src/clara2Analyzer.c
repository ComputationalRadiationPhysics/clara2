/**
 * Copyright 2014 Alexander Debus, Lucas Clemente (flint)
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

///\file clara2Analyzer.c included into GPT.
///A GPT custom element, that launches, initializes, calls and finalizes Clara2 within GPT.

#include "clara2Analyzer.h"

///The time of the previous timestep
double t_old;

///Kill Clara2
void exitClara2(void* p) {
	outputData();
}

///Update the 'pars_old' array
void copyArray()
{
	int i;
	for (i = 0; i < numpar; i++)
	{
		old_par* p = &(pars_old[pars[i].ID]);
		p->GBx = pars[i].GBr[0];
		p->GBy = pars[i].GBr[1];
		p->GBz = pars[i].GBr[2];
		p->G = pars[i].G;
		p->stored = 1;
	}
	//New step
	step++;
}

///Called after every succesful timestep
int calcClara2(double t, double* dt, double* x, void* info);

/* Info structure containing all relevant parameters for this element */
struct clara2Analyzer_info
{
} ;

///The init function of flint.
void clara2Analyzer_init(gptinit *init) {
  
  struct clara2Analyzer_info *info ;

  /* Read Element Coordinate System (ECS) from parameter list */
  gptbuildECS( init ) ;

  /* Print usage line when the number of parameters is incorrect */
  if( gptgetargnum(init)!=0 )
    gpterror( "Syntax: %s(ECS)\n", gptgetname(init) ) ;
  
  /* Allocate memory for info structure */
  info = (struct clara2Analyzer_info*)gptmalloc( sizeof(struct clara2Analyzer_info) ) ;
  
  printf("GPT: clara2Analyzer_init.");
  
  //Allocate and set the 'pars_old' array
  ///\todo Not very nice, but neccesary, as numpars can change during simulation
  pars_old = (old_par*)gptmalloc(sizeof(old_par) * NUMPARS_SUPPORTED);
  copyArray();
  int i;
  for (i = 0; i < NUMPARS_SUPPORTED; i++) {
    pars_old[i].stored = 0;
  }
  
  //Initialize step number and time
  step = 0;
  t_old = 0; // To do: Initialize it with tstart
  
  //Of course...
  pars_clara2 = 0;
  
  //Start Clara2
  if (!clara2_init()) gpterror("Clara2 could not init. Aborting...\n");
  
  //Exit-function
  gptaddmainfunction(GPTMAINFNC_EXT, &exitClara2, 0);
  //Main-function
  odeaddoutfunction(ODEFNC_USR, &calcClara2, 0);
}

///Main simulation routine in GPT
int calcClara2(double t, double* dt, double* x, void* info)
{
	///\todo Hopefully temporally... See NUMPARS_SUPPORTED.
	if (numpar >= NUMPARS_SUPPORTED)
	{
		gpterror("The number of particles (=%d) has exceeded the maximum (=%d). Large numbers are not implemented yet!\n", numpar, NUMPARS_SUPPORTED);
	}
  
	//First step? Do nothing.
	if (step == 0)
	{
		printf("First step...\n");
		copyArray();
		return 0;
	}
  	
	if ( !clara2_calculate(numpar, t ,t - t_old) )
	{
		gpterror("Clara2 exited with error. It should have produced output...\n");
	}
	
	//Print a status message once in a while
	if ((step % OUTPUT_INTERVAL) == 0)
	{
		printf("Processing... step=%d, t=%e, dt=%e, clara2=%.2f, alive=%d\n", step, t, *dt,
			(double)pars_clara2 / (double) alive * 100.0, alive);
	}

	//Update old-pars
	copyArray();
	//Update time information
	t_old = t;
  	return 0;
}
