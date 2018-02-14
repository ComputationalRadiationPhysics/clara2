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




template <typename X, typename Y>
void interpolation(const X* x_old,
                   const Y* y_old,
                   const unsigned N_old,
                   const X* x_new,
                   Y* y_new,
                   const unsigned N_new)
{
  for(unsigned i_new = 0, i_old = 0; i_new< N_new; ++i_new)
  {
    while(x_old[i_old] < x_new[i_new] && i_old < N_old)
      ++i_old;
    if(i_old < N_old && i_old > 0)
    {
      y_new[i_new] = (Y)((y_old[i_old] - y_old[i_old -1])/(x_old[i_old] - x_old[i_old -1])
                         * (x_new[i_new] - x_old[i_old -1]) + y_old[i_old -1]);
    }
    else
    {
      y_new[i_new] = 0;
    }
  }
}


template <typename X, typename Y>
void interpolation_on(const X* x_old,
                      const Y* y_old,
                      const unsigned N_old,
                      const X* x_new,
                      Y* y_new,
                      const unsigned N_new)
{
  for(unsigned i_new = 0, i_old = 0; i_new< N_new; ++i_new)
  {
    while(x_old[i_old] < x_new[i_new] && i_old < N_old)
      ++i_old;
    if(i_old < N_old && i_old > 0)
    {
      y_new[i_new] += (Y)((y_old[i_old] - y_old[i_old -1])/(x_old[i_old] - x_old[i_old -1])
                          * (x_new[i_new] - x_old[i_old -1]) + y_old[i_old -1]);
    }
    else
    {
      y_new[i_new] += (Y)0;  // could be skipped
    }
  }
}


void interpolation_on(Detector_fft* fft,
                      const double* x_new,
                      double* y_new,
                      const unsigned N_new)
{
  unsigned N_old = fft->half_frequency();

  for(unsigned i_new = 0, i_old = 0; i_new< N_new; ++i_new)
  {
    while(fft->get_spectrum(i_old, 0) < x_new[i_new] && i_old < N_old)
      ++i_old;
    if(i_old < N_old && i_old > 0)
    {
      y_new[i_new] += ((fft->get_spectrum(i_old, 1) - fft->get_spectrum(i_old-1, 1))/(fft->get_spectrum(i_old, 0) - fft->get_spectrum(i_old-1, 0))
                       * (x_new[i_new] - fft->get_spectrum(i_old-1, 0)) + fft->get_spectrum(i_old-1, 1));
    }
    else
    {
      y_new[i_new] += 0;  // could be skipped
    }
  }
}


void interpolation_int(Detector_fft* fft,
                       const double* x_new,
                       double* y_new,
                       const unsigned N_new)
{
  unsigned N_old = fft->half_frequency();
  unsigned total_counter =0;
  unsigned j=0;

  //in between:
  for(unsigned  i=0; i<N_new; ++i)
  {
    double x_high;
    double x_low;
    if(i==0)
    {
      x_high = 0.5 * (x_new[i] + x_new[i+1]);
      x_low  = x_new[i];
    }
    else if(i==N_new-1)
    {
      x_high = x_new[i];
      x_low  = 0.5 * (x_new[i-1] +  x_new[i]);
    }
    else
    {
      x_high = 0.5 * (x_new[i] + x_new[i+1]);
      x_low  = 0.5 * (x_new[i-1] +  x_new[i]);
    }


    unsigned counter=0;
    double sum = 0.0;


    while(fft->get_spectrum(j, 0) >= x_low && fft->get_spectrum(j, 0) <= x_high && j<N_old)
    {
      counter++;
      sum += fft->get_spectrum(j, 1);
      ++j;
    }

    if (counter == 0) // this avoids a NaN if the bin does not contain any data
    {
      y_new[i] += 0.0;
    }
    else
    {
      y_new[i] += sum/counter;
    }

    total_counter += counter;

  }
}
