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

#include "physics_units.hpp"
#include "utilities.hpp"


// Constructors and Destructors:

/** constructor for a spectral detector
 *
 * @param n_unit      = unit vector pointing in observation direction
 * @param spec_length = number of frequencies to compute
 * @param omega_max   = maximum frequency -> frequency range [0, omega_max]
 */
Detector_dft::Detector_dft(const R_vec n_unit,
                           const unsigned spec_length,
                           const double omega_max)
  : n_unit(n_unit.unit_vec()),
    spec_length(spec_length),
    spectrum(0),
    spectrum_mag(0),
    frequency(0)
{
  spectrum = new Vector<std::complex<double>, 3>[spec_length]; /* complex amplitude for each omega */
  spectrum_mag = new double[spec_length]; /* absolute square of spectrum */
  frequency = new double[spec_length]; /* omega values */

  set_frequency(omega_max); /* compute frequencies */
}

/** constructor for a spectral detector
 *
 * @param n_unit      = unit vector pointing in observation direction
 * @param spec_length = number of frequenies to compute
 * @param omega       = pointer to frequencies to calculate
 */
Detector_dft::Detector_dft(const R_vec n_unit,
                           const unsigned spec_length,
                           const double* omega)
  : n_unit(n_unit.unit_vec()),
    spec_length(spec_length),
    spectrum(0), spectrum_mag(0),
    frequency(0)
{
  spectrum = new Vector<std::complex<double>, 3>[spec_length]; /* complex amplitude for each omega */
  spectrum_mag = new double[spec_length]; /* absolute square of spectrum */
  frequency = new double[spec_length]; /* omega values */

  set_frequency(omega); /* copy frequency values */
}

/** destructor */
Detector_dft::~Detector_dft()
{
  delete[] spectrum;
  delete[] spectrum_mag;
  delete[] frequency;
}


// add to spectrum methods:

/** method to add an radiation amplitude for one time step
 *
 * @param r_0 = position
 * @param beta_0 = speed / speed_of_light
 * @param dot_beta_0 = time derivative of beta_0
 * @param t_part_0 = time
 * @param delta_t = time step width
 */
void Detector_dft::add_to_spectrum(const R_vec r,
                                   const R_vec beta,
                                   const R_vec dot_beta,
                                   const double t_part,
                                   const double delta_t)
{
  /* calculate (real) vector part of radiation amplitude */
  const R_vec fou1 = (n_unit%((n_unit-beta)%dot_beta))
                      / util::square(1. - beta * n_unit);

  /* make vector part complex for calculation in complex space */
  const Vector<std::complex<double>, 3> fou1_complex = fou1.make_complex();

  /* for each frequency: */
  for(unsigned a = 0; a < spec_length; ++a)
  {
    const double t_ret = t_part - (n_unit*r)/phy::c; /* retarded time */
    /* compute complex phase for this retarded time */
    std::complex<double> fou2 = std::polar(1.0, frequency[a]*t_ret);
    /* for each dimension: */
    for (unsigned i=0; i<3; ++i)
    {
      /* this is the (local) Nyquist limiter - only add amplitude if
       * frequency is below (local) Nyquist frequency */
      if(frequency[a] < 0.75 * M_PI / (delta_t*(1.0 - beta * n_unit)) )
        (spectrum[a])[i] += fou1_complex[i]*fou2*delta_t;
    }

  }
}

/* ISSUE # 74 - interface too crowded */
/** method to add an radiation amplitude for one time step
 *
 * @param r_0 = position
 * @param p_0 = momentum
 * @param dot_p_0 = time derivative of p_0
 * @param beta_0 = speed / speed_of_light
 * @param gamma_0 = relativistic gamma factor
 * @param dot_gamma_0 = time derivative of gamma_0
 * @param t_part_0 = time
 * @param delta_t = time step width
 */
void Detector_dft::add_to_spectrum(const R_vec r,
                                   const R_vec p,
                                   const R_vec dot_p,
                                   const R_vec beta,
                                   const double gamma,
                                   const double dot_gamma,
                                   const double t_part,
                                   const double delta_t)
{
  /* normalize momentum */
  const R_vec     p_wave = p/(phy::m_e*phy::c);
  const R_vec dot_p_wave = dot_p/(phy::m_e*phy::c);

  /* calculate (real) vector part of radiation amplitude */
  const R_vec fou1 = (n_unit%((gamma*n_unit - p_wave)%(dot_p_wave - beta*dot_gamma)))
                     / util::square(gamma - p_wave*n_unit);

  /* make vector part complex for calculation in complex space */
  const Vector<std::complex<double>, 3> fou1_complex = fou1.make_complex();

  /* for each frequency: */
  for(unsigned a = 0; a < spec_length; ++a)
  {
    const double t_ret = t_part - (n_unit*r)/phy::c; /* retarded time */
    /* compute complex phase for this retarded time */
    std::complex<double> fou2 = std::polar(1.0, frequency[a]*t_ret);
    /* for each dimension: */
    for (unsigned i=0; i<3; ++i)
    {
      /* this is the (local) Nyquist limiter - only add amplitude if
       * frequency is below (local) Nyquist frequency */
      if(frequency[a] < 0.75 * M_PI / (delta_t*(1.0 - beta * n_unit)) )
        (spectrum[a])[i] += fou1_complex[i]*fou2*delta_t;
    }

  }
}


// calculate spectrum method:

/** calculate the spectrum from all added amplitudes */
void Detector_dft::calc_spectrum()
{
  /* SI unit factor for calculating radiation energy per frequency and solid angle */
  const double factor = util::square(phy::q)
                        / (16.*util::cube(M_PI)*phy::epsilon_0*phy::c);

  /* calculate the absolute value for each frequency */
  for(unsigned a = 0; a < spec_length; ++a)
  {
    spectrum_mag[a] = factor * util::square<R_vec, double>(spectrum[a].abs());
  }
}


// Getter:

/* ISSUE #75 - separate get spectrum and get frequency */
/** get spectrum or frequency
  *
  * @param a = id between [0, number of frequencies]
  * @param b = int return 0->frequency, 1->spectra in Js
  */
double Detector_dft::get_spectrum(unsigned a,
                                  unsigned b)
{
  assert(a < spec_length);
  switch(b)
  {
    case 0: /* return frequency */
      return frequency[a];
      break;
    case 1: /* return spectral value at that frequency */
      return spectrum_mag[a];
      break;
    default: /* index out of range */
      std::cout << "Wrong access to spectrum (dft, beta)!" << std::endl;
      std::cout << b << " is larger than 2." << std::endl;
      assert(false);
      break;
  }
}

/** return total energy (in calculated spectral range) */
double Detector_dft::energy()
{
  double result = 0.; /* total energy */
  /* sum over all frequencies */
  for (unsigned i = 0; i < (spec_length); ++i)
    result += spectrum_mag[i];

  /* multiply by frequency bin width */
  /* ISSUE #76 - remove this arbitrary delta_omega here */
  result *= (frequency[7] - frequency[6]);
  return result;
}


// set frequency methods:

/** set frequency based on maximum frequency
  * to cover frequency range [0, omega_max]
  *
  * @ param omega_max = maximum frequency
  */
inline void Detector_dft::set_frequency(const double omega_max)
{
  for (unsigned i=0; i<spec_length; ++i)
    frequency[i] = (double(i))/(double(spec_length)) * omega_max;
}

/** set frequency using pointer to frequency values
  *
  * @ param omega = double pointer to frequency values
  */
inline void Detector_dft::set_frequency(const double* omega)
{
  /* copy frequency values */
  for (unsigned i=0; i<spec_length; ++i)
    frequency[i] =  omega[i];
}
