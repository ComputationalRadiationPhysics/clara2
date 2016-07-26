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



#include <cmath>
#include "utilities.hpp"
#include "physics_units.hpp"

#ifndef ANALYTICAL_SOLUTION_RPAUSCH
#define ANALYTICAL_SOLUTION_RPAUSCH



namespace analy
{


const int N_per = 225;
const double lambda_u = 0.4E-6;
const double Bmax = 1.0;
const unsigned n_mode = 1;

double K;
double F_approx;


using namespace std;

double L(double x)
{
  return (util::square(sin(M_PI*x))) / (util::square(N_per*sin(M_PI*x/N_per)));
	  
}

double L2(double x)
{
  return (util::square(sin(M_PI*x))) / (util::square(M_PI*x));
	  
}



double Lambda(double gamma)
{
  return lambda_u/(2*util::square(gamma)) * (1. + util::square(K)/2.);
}

double Spectrum(double omega, double gamma)
{
  return util::square(phy::q*N_per*gamma)/
    (4.*M_PI*phy::epsilon_0*phy::c) * 
    L(N_per*(omega-2.*M_PI*phy::c/Lambda(gamma))/(2.*M_PI*phy::c / 
						  Lambda(gamma))  )
    * F_approx;
}

double Spectrum2(double omega, double gamma)
{
  return util::square(phy::q*N_per*gamma)/
    (4.*M_PI*phy::epsilon_0*phy::c) * 
    L2(N_per*(omega-2.*M_PI*phy::c/Lambda(gamma))/(2.*M_PI*phy::c / 
						  Lambda(gamma))  )
    * F_approx;
}



void calc_K(double gamma)
{
  K = 0.0569; //phy::q * Bmax * lambda_u / (2. * M_PI * phy::m_e * phy::c);
  F_approx = util::square(K * n_mode) / util::square(1. + util::square(K)/2.) ;

  std::cout << " K= " << K << "\t F_approx = " << F_approx << std::endl;
}


}

#endif
