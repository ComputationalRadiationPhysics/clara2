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



#include "detector_dft.hpp"


// Constructors and Destructors:

Detector_dft::Detector_dft(R_vec n_unit, const unsigned spek_length,
			   const double omega_max)
  : n_unit(n_unit.unit_vec()),  spek_length(spek_length),
    spektrum(0), spektrum_mag(0), frequency(0)
{
  spektrum = new Vector<std::complex<double>, 3>[spek_length];
  spektrum_mag = new double[spek_length];
  frequency = new double[spek_length];

  set_frequency(omega_max); 
}

Detector_dft::Detector_dft(R_vec n_unit,  const unsigned spek_length,
			   const double* omega)
  : n_unit(n_unit.unit_vec()),  spek_length(spek_length),
    spektrum(0), spektrum_mag(0), frequency(0)
{
  spektrum = new Vector<std::complex<double>, 3>[spek_length];
  spektrum_mag = new double[spek_length];
  frequency = new double[spek_length];

  set_frequency(omega); 
}


Detector_dft::~Detector_dft()
{
  delete[] spektrum;
  delete[] spektrum_mag;
  delete[] frequency;
}


// add to spectrum methods:

void Detector_dft::add_to_spectrum(const R_vec r, 
				   const R_vec beta,
				   const R_vec dot_beta,
                   const double t_part,
                   const double delta_t)
{
  const R_vec fou1a = (n_unit%((n_unit-beta)%dot_beta)) 	
                        / util::square(1. - beta * n_unit);

 
  const Vector<std::complex<double>, 3> fou1b = fou1a.make_complex();

  for(unsigned a = 0; a < spek_length; ++a)
    {
      std::complex<double> fou2 = std::polar(1.0, frequency[a]*
		                  (t_part - (n_unit*r)/phy::c));
      for (unsigned i=0; i<3; ++i)
	{
        // this is the (local) Nyquist limiter 
        if(frequency[a] < 0.75 * M_PI / (delta_t*(1.0 - beta * n_unit)) )
	  (spektrum[a])[i] += fou1b[i]*fou2*delta_t;
	}

    }
}

void Detector_dft::add_to_spectrum(const R_vec r, 
				   const R_vec p, 
				   const R_vec dot_p,
				   const R_vec beta,
				   const double gamma,
				   const double dot_gamma,
                   const double t_part,
                   const double delta_t)
{		
  const R_vec     p_wave = p/(phy::m_e*phy::c);
  const R_vec dot_p_wave = dot_p/(phy::m_e*phy::c);

  const R_vec fou1 = (n_unit%((gamma*n_unit - p_wave)%(dot_p_wave - beta*dot_gamma)))
	               / util::square(gamma - p_wave*n_unit);

  const Vector<std::complex<double>, 3> fou1_complex = fou1.make_complex();

  for(unsigned a = 0; a < spek_length; ++a)
    {
      std::complex<double> fou2 = std::polar(1.0, frequency[a]*
				  (t_part - (n_unit*r)/phy::c));
      for (unsigned i=0; i<3; ++i)
	{
        // this is the (local) Nyquist limiter 
        if(frequency[a] < 0.75 * M_PI / (delta_t*(1.0 - beta * n_unit)) )
	  (spektrum[a])[i] += fou1_complex[i]*fou2*delta_t;
	}

    }
}


// calculate spectrum method:

void Detector_dft::calc_spectrum()
{
  const double factor = util::square(phy::q)/
                (16.*util::cube(M_PI)*phy::epsilon_0*phy::c);

  for(unsigned a = 0; a < spek_length; ++a)
    {
      spektrum_mag[a] = factor *
	                util::square<R_vec, double>(spektrum[a].abs());
    }
}


// Getter:

double Detector_dft::get_spectrum(unsigned a, unsigned b)
{
  assert(a < spek_length);
  switch(b)
    {
    case 0:
      return frequency[a];
      break;
    case 1:
      return spektrum_mag[a];
      break;
    default:
      std::cout << "Wrong access to spectrum (dft, beta)!" << std::endl;
      std::cout << b << " is larger than 2." << std::endl;
      assert(false);
      break;
    }
}

double Detector_dft::energy()
{
  double result = 0.;
  for (unsigned i = 0; i < (spek_length); ++i)
    result += spektrum_mag[i];
  result *= (frequency[7] - frequency[6]);
  return result;
}


// set frequency methods:

inline void Detector_dft::set_frequency(const double omega_max)  
{
  for (unsigned i=0; i<spek_length; ++i) 
    frequency[i] = (double(i))/(double(spek_length)) * omega_max;
}

inline void Detector_dft::set_frequency(const double* omega)  
{
  for (unsigned i=0; i<spek_length; ++i) 
    frequency[i] =  omega[i];
}



