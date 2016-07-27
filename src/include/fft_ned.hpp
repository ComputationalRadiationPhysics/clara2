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





#include <iostream>
#include <fstream>
#include <cmath>
#include <time.h>
#include <stdio.h>
#include <complex>



#ifndef NED_FFT_RPAUSCH
#define NED_FFT_RPAUSCH

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
ned_FFT(unsigned N, A x_original[], T y_original[])
  : N_data(N), x_equi(0), y_equi(0), data_complex(0), spektrum(0), omega(0)
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
inline void swap(T& a, T& b)
{
  T temp = a;
  a = b;
  b = temp;
}

void fft(T* data, unsigned long nn)
{
    unsigned long n, mmax, m, j, istep, i;
    A wtemp, wr, wpr, wpi, wi, theta;  
    // d     d    d     d   d   d
    T tempr, tempi;
    // data  data
 
    // reverse-binary re-indexing
    n = nn<<1;
    j=1;
    for (i=1; i<n; i+=2) {
        if (j>i) {
            swap(data[j-1], data[i-1]);
            swap(data[j], data[i]);
        }
        m = nn;
        while (m>=2 && j>m) {
            j -= m;
            m >>= 1;
        }
        j += m;
    };
 
    // here begins the Danielson-Lanczos section
    mmax=2;
    while (n>mmax) {
        istep = mmax<<1;
        theta = -(2*M_PI/mmax);
        wtemp = std::sin(0.5*theta);
        wpr = -2.0*wtemp*wtemp;
        wpi = std::sin(theta);
        wr = 1.0;
        wi = 0.0;
        for (m=1; m < mmax; m += 2) {
            for (i=m; i <= n; i += istep) {
                j=i+mmax;
                tempr = wr*data[j-1] - wi*data[j];
                tempi = wr * data[j] + wi*data[j-1];
 
                data[j-1] = data[i-1] - tempr;
                data[j] = data[i] - tempi;
                data[i-1] += tempr;
                data[i] += tempi;
            }
            wtemp=wr;
            wr += wr*wpr - wi*wpi;
            wi += wi*wpr + wtemp*wpi;
        }
        mmax=istep;
    }
}

A interpolation_equi(A x_0[], T y_0[], unsigned N_0, 
		     A x_1[], T y_1[], unsigned N_1)
{
  for (unsigned i=1; i < N_0; ++i)
    if(x_0[i-1] > x_0[i]) 
      {std::cerr << "error 01: interpolation inverted (fft_ned.hpp) " << i << " --> " << x_0[i-1]
		 << " <=! " << x_0[i] << "\n";}

  const A min = x_0[0];
  const A max = x_0[N_0-1];

  // creating equidistant x_values
  for (unsigned i=0; i < N_1; ++i)
    x_1[i] = min + (max-min)/(N_1) *i ;


  // calculating y_values
  unsigned j=0;
  for (unsigned i=0; i<N_1; ++i)
    {
      for(; !(x_0[j] <= x_1[i] && x_1[i] < x_0[j+1]); ++j) {}

      if (!(x_0[j] <= x_1[i] && x_1[i] < x_0[j+1]))
	std::cerr << "error 02: (fft_ned.hpp)" << std::endl;

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




public:
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




#endif

