/**
 * Copyright 2014-2017 Richard Pausch
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

#include <fftw3.h>
#pragma once



inline unsigned power_of_two(unsigned N)
{
  unsigned exponent=1;
  for(; N > (1u<<exponent); ++exponent) {}
  return (1u<<exponent);
}


// non equal distant FFT
template< typename A, typename T >  // A...time,   T...data                                                                                                                           
class ned_FFT
{

public:
  // constructor
  ned_FFT(unsigned N,
          A x_original[],
          T y_original[])
    : N_data(N),
      x_equi(0),
      y_equi(0),
      data_complex(0),
      spektrum(0),
      omega(0)
  {
    x_equi = new A[N];
    y_equi = new T[N];

    delta_t = interpolation_equi(x_original, y_original,  N_data,
                                 x_equi, y_equi, N_data);

    unsigned exponent=1;
    for(; N_data > (1u<<exponent); ++exponent) {}
    N_bin = 1u<<exponent;

    data_complex = new T[N_bin<<1];
    for(unsigned i=0; i<N_data; ++i)
    {
      data_complex[2*i] = y_equi[i];
      data_complex[2*i+1] = T(0.);
    }

    for(unsigned i=2*N; i<(N_bin<<1); ++i)
      data_complex[i] = T(0.);

    fft(data_complex, N_data);

    omega_calc();
    spektrum_calc();
  }

  // destructor
  ~ned_FFT()
  {
    delete[] x_equi;
    delete[] y_equi;
    delete[] data_complex;

    delete[] spektrum;
    delete[] omega;
  }

private:

  void fft(T* data,
           unsigned long N)
  {
    // this uses the awesome fftw3 library
    // see http://www.fftw.org/

    // transfer data to fftw3 own data structure
    // TODO we could avoid using this data transfer
    // by using fftw data types right away
    // for now not, to keep fft independent to allow later
    // use of liFFT
    // https://github.com/ComputationalRadiationPhysics/liFFT
    fftw_complex *input, *output;
    input = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    output = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

    // default fft complex to complex
    // TODO since input is real, there is a faster real implementation
    // without Nyquist reflections
    fftw_plan plan =  fftw_plan_dft_1d(N,
                                       input,
                                       output,
                                       FFTW_FORWARD,
                                       FFTW_ESTIMATE);
    
    // here we have to assume that the data is a vector type of 3 dimensions
    // TODO this is not necessarily the case
    // and this breaks the template structure of the rest of the code
    for(unsigned int index_vec = 0; index_vec < 3; index_vec++)
      {
        for (unsigned int i=0; i<N; i++)
          {
            input[i][0] = data[2*i][index_vec]; // real input (signal)
            input[i][1] = data[2*i+1][index_vec]; // imaginary input (signal)
          }
    
        fftw_execute(plan); // run FFT

        // copy data back into original data array
        for(unsigned int i=0; i<N; i++)
          {
            data[2*i][index_vec] = output[i][0]; // real output (spectrum)
            data[2*i+1][index_vec] = output[i][1]; // imaginary output (spectrum)
          }
      }
    // free memory for fftw in-between data
    fftw_destroy_plan(plan); // could be freed before data transfer
    fftw_free(input); fftw_free(output);
  }

  A interpolation_equi(A x_0[],
                       T y_0[],
                       unsigned N_0,
                       A x_1[],
                       T y_1[],
                       unsigned N_1)
  {
    for (unsigned i=1; i < N_0; ++i)
    {
      if(x_0[i-1] > x_0[i])
      {
        std::cerr << "error 01: interpolation inverted (ned_fft.hpp) "
                  << i << " --> " << x_0[i-1]
                  << " <=! " << x_0[i] << "\n";
      }
    }

    const A min = x_0[0];
    const A max = x_0[N_0-1];

    // creating equidistant x_values
    for (unsigned i=0; i < N_1; ++i)
      x_1[i] = min + (max-min)/(N_1) * i;

    // calculating y_values
    unsigned j=0;
    for (unsigned i=0; i<N_1; ++i)
    {
      for(; !(x_0[j] <= x_1[i] && x_1[i] < x_0[j+1]); ++j) {}

      if (!(x_0[j] <= x_1[i] && x_1[i] < x_0[j+1]))
        std::cerr << "error 02: (ned_fft.hpp)" << std::endl;

      y_1[i] = y_0[j] + (y_0[j+1] - y_0[j])*((x_1[i]-x_0[j])/(x_0[j+1]-x_0[j]));
    }
    return (max-min)/N_1;
  }


// calculate angular frequency
  void omega_calc()
  {
    omega = new A[N_bin];
    for (unsigned i=0; i<N_bin; ++i)
      omega[i] = (2*M_PI*i)/(N_bin*delta_t);
  }

  // calculate spectrum (T-->A)
  void spektrum_calc()
  {
    spektrum = new A[N_bin];
    for (unsigned i=0; i<N_bin; ++i)
      spektrum[i] = std::sqrt(data_complex[2*i]*data_complex[2*i] +
                              data_complex[2*i+1]*data_complex[2*i+1]);
  }


public:
  unsigned N_data;
  unsigned N_bin;
  A delta_t;
  A* x_equi;
  T* y_equi;
  T* data_complex;
  A* spektrum;
  A* omega;

};



// usage:
//  for (unsigned i = 0; i< N_data; ++i)
//    {
//      x_data[i] = random_double(i, 0.0)*0.2 - 5.003;
//      y_data[i] = fkt(x_data[i]);
//    }
//
//  ned_FFT<double> spektrum(N_data, x_data, y_data);
