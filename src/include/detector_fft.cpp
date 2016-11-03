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


// Constructor and Destructor:
Detector_fft::Detector_fft(R_vec n_unit,
                           unsigned N_data)
  : n_unit(n_unit.unit_vec()),
    N_data(N_data),
    counter(0),
    time(0),
    data(0),
    spektrum_mag(0),
    frequency(0)
{
  spek_length = power_of_two(N_data)*fft_length_factor;

  //std::cout << "spek_length : " << spek_length << std::endl;
  time = new double[spek_length];
  data = new R_vec[spek_length];
  spektrum_mag = new double[spek_length];
  frequency = new double[spek_length];
}

Detector_fft::~Detector_fft()
{
  delete[] time;
  delete[] data;
  delete[] spektrum_mag;
  delete[] frequency;
}


// add to spectrum methods:

void Detector_fft::add_to_spectrum(const R_vec r,
                                   const R_vec beta,
                                   const R_vec dot_beta,
                                   const double t_part,
                                   const double delta_t)
{
  const R_vec fou1 = (n_unit%((n_unit-beta)%dot_beta))
                     / util::cube(1. - beta * n_unit);

  assert(counter < spek_length); // --> still necessary? probably not

  time[counter] = t_part - (n_unit*r)/phy::c;
  data[counter] = fou1;

  ++counter;
}


void Detector_fft::add_to_spectrum(const R_vec r,
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

  const R_vec fou = (n_unit%((gamma*n_unit - p_wave)%(dot_p_wave
                      - beta*dot_gamma))) / util::cube(gamma - p_wave*n_unit)
                     * gamma;

  assert(counter < spek_length); // --> still necessary? probably not

  time[counter] = t_part - (n_unit*r)/phy::c;
  data[counter] = fou;

  ++counter;
}


// calculate spectrum method:

void Detector_fft::calc_spectrum()
{
  for (unsigned i= counter; i<spek_length; ++i)
    time[i] = time[counter-1];

  const double factor = util::square(phy::q)
                        / (16.*util::cube(M_PI)*phy::epsilon_0*phy::c);

  ned_FFT<double, R_vec> analyse(spek_length, time, data);

  for(unsigned i=0; i<spek_length; ++i)
  {
    frequency[i] = analyse.omega[i];
    spektrum_mag[i] = factor * util::square(analyse.spektrum[i] * analyse.delta_t);
  }
}


double Detector_fft::get_spectrum(unsigned a, unsigned b)
{
  assert(a < (spek_length)   );
  switch(b){
  case 0:
    return frequency[a];
    break; //necessary?
  case 1:
    return spektrum_mag[a];
    break;
  default:
    std::cerr << "Wrong access to spectrum (fft, momentum)!" << std::endl;
    std::cerr << b << " is larger than 1." << std::endl;
    assert(false);
    break;
  }
}


double Detector_fft::energy()
{
  double result = 0;
  for (unsigned i = 0; i < half_frequency(); ++i)
    result += spektrum_mag[i];

  result *= (frequency[7] - frequency[6]);
  return result;
}


unsigned Detector_fft::half_frequency()
{
  return spek_length>>1;
}
