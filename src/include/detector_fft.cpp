/**
 * Copyright 2014-2016 Richard Pausch, Alexander Koehler
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


#include "detector_fft.hpp"

#include "../settings.hpp"
#include "physics_units.hpp"
#include "utilities.hpp"
#include "ned_fft.hpp"


// Constructor and Destructor:

/** constructor for a point-like detector using a Fast Fourier Transform
  * to calculate the radiation spectra
  *
  *  @param n_unit   = unit vector in direction of energy deposition
  *  @param N_data   = number of time step in trace
  */
Detector_fft::Detector_fft(const R_vec n_unit,
                           const unsigned N_data)
  : n_unit(n_unit.unit_vec()),
    N_data(N_data),
    counter(0),
    time(0),
    data(0),
    spectrum_mag(0),
    frequency(0)
{
  spec_length = power_of_two(N_data) * param::fft_length_factor;

  time = new double[spec_length]; /* retarded time */
  data = new R_vec[spec_length]; /* real vector amplitude */
  spectrum_mag = new double[spec_length]; /* absolute square of spectrum */
  frequency = new double[spec_length]; /* omega values */
}

/** destructor */
Detector_fft::~Detector_fft()
{
  delete[] time;
  delete[] data;
  delete[] spectrum_mag;
  delete[] frequency;
}


// add to spectrum methods:

/** method to add an radiation amplitude for one time step
  *
  * @param r = position
  * @param beta = speed / speed_of_light
  * @param dot_beta = time derivative of beta
  * @param t_part = time
  * @param delta_t = time step width
  */
void Detector_fft::add_to_spectrum(const R_vec r,
                                   const R_vec beta,
                                   const R_vec dot_beta,
                                   const double t_part,
                                   const double delta_t)
{
  /* calculate (real) vector part of radiation amplitude */
  const R_vec fou1 = (n_unit%((n_unit-beta)%dot_beta))
                     / util::cube(1. - beta * n_unit);

  /* abort if trace is longer than spectra */
  assert(counter < spec_length); // --> still necessary? probably not

  /* retarded time */
  time[counter] = t_part - (n_unit*r)/phy::c;
  /* real vector amplitude */
  data[counter] = fou1;

  ++counter; /* count processed time steps */
}


/* ISSUE # 74 - interface too crowded */
/** method to add an radiation amplitude for one time step
  *
  * @param r = position
  * @param p = momentum
  * @param dot_p = time derivative of p
  * @param beta = speed / speed_of_light
  * @param gamma = relativistic gamma factor
  * @param dot_gamma = time derivative of gamma
  * @param t_part = time
  * @param delta_t = time step width
  */
void Detector_fft::add_to_spectrum(const R_vec r,
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
  const R_vec fou = (n_unit%((gamma*n_unit - p_wave)%(dot_p_wave
                      - beta*dot_gamma))) / util::cube(gamma - p_wave*n_unit)
                     * gamma;

  /* abort if trace is longer than spectra */
  assert(counter < spec_length); // --> still necessary? probably not

  /* retarded time */
  time[counter] = t_part - (n_unit*r)/phy::c;
  /* real vector amplitude */
  data[counter] = fou;

  ++counter; /* count processed time steps */
}


// calculate spectrum method:

/** calculate the spectrum from all added amplitudes */
void Detector_fft::calc_spectrum()
{
  /* fill the last entries between N_data and spec_length
   * with the last time = results in zero-padding */
  for (unsigned i= counter; i<spec_length; ++i)
    time[i] = time[counter-1];

  /* SI unit factor for calculating radiation energy per frequency and solid angle */
  const double factor = util::square(phy::q)
                        / (16.*util::cube(M_PI)*phy::epsilon_0*phy::c);

  /* compute spectra using a not-equidistant FFT method */
  ned_FFT<double, R_vec> analyse(spec_length, time, data);

  /* calculate the absolute value for each frequency */
  for(unsigned i=0; i<spec_length; ++i)
  {
    frequency[i] = analyse.omega[i];
    spectrum_mag[i] = factor * util::square(analyse.spektrum[i] * analyse.delta_t);
  }
}


/* ISSUE #75 - separate get spectrum and get frequency */
/** get spectrum or frequency
  *
  * @param a = id between [0, number of frequencies]
  * @param b = int return 0->frequency, 1->spectra in Js
  */
double Detector_fft::get_spectrum(unsigned a, unsigned b)
{
  assert(a < (spec_length)   );
  switch(b){
  case 0: /* return frequency */
    return frequency[a];
    break; //necessary?
  case 1: /* return spectral value at that frequency */
    return spectrum_mag[a];
    break;
  default: /* index out of range */
    std::cerr << "Wrong access to spectrum (fft, momentum)!" << std::endl;
    std::cerr << b << " is larger than 1." << std::endl;
    assert(false);
    break;
  }
}

/** return total energy (in calculated spectral range) */
double Detector_fft::energy()
{
  double result = 0; /* total energy */
  /* sum over all frequencies */
  for (unsigned i = 0; i < half_frequency(); ++i)
    result += spectrum_mag[i];

  /* multiply by frequency bin width */
  /* ISSUE #76 - remove this arbitrary delta_omega here */
  result *= (frequency[7] - frequency[6]);
  return result;
}


/** return the number of frequencies till the Nyquist frequency
  * (number of half of all frequency bins)
  */
unsigned Detector_fft::half_frequency()
{
  return spec_length>>1; /* = spec_length / 2 */
}
