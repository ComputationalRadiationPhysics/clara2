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



#include "vector.hpp"


// constructor, destructor:

//basic constructor
template< typename T, unsigned N>
inline Vector<T, N>::Vector()     
{
  for (unsigned i=0; i<N; ++i) {
    data[i] = T(0.0);
  }
}

// constructor for 3D
template< typename T, unsigned N>
inline Vector<T, N>::Vector(double x, double y, double z)
{
  assert( N==3 ); // better via template? one more dummy function
  data[0]=x;
  data[1]=y;
  data[2]=z;
}

       


////////////////////////////////////////////////
// Getters and Setters:

// access data
template< typename T, unsigned N>
inline T & Vector<T, N>::operator[](unsigned i) 
{
  assert(i<N);
  return data[i];
}

// access data
template< typename T, unsigned N> 
const inline T & Vector<T, N>::operator[](unsigned i) const 
{
  assert(i<N);
  return data[i];
}



/////////////////////////////////////////////////////   
// Calculations:
    
/// addition
template< typename T, unsigned N> 
inline Vector<T, N> Vector<T, N>::operator+ (const  Vector& v) const  
{
  Vector back;
  for (unsigned i=0; i<N; ++i) 
    {
      back[i] = data[i] + v.data[i];
    }
  return back;
}

/// subtraction
template< typename T, unsigned N> 
inline Vector<T, N> Vector<T, N>::operator- (const Vector& v) const  
{
  Vector back;
  for (unsigned i=0; i<N; ++i) 
    {
      back[i] = data[i] - v.data[i];
    }
  return back;
}
    
/// magnitude
template< typename T, unsigned N> 
inline T Vector<T, N>::mag() const  
{ 
  T sqr_magnitude=0;
  for (unsigned i=0; i<N; ++i) 
    {
      sqr_magnitude += data[i] * data[i];
    }
  return std::sqrt(sqr_magnitude); 
}

/// scalar product
template< typename T, unsigned N> 
inline T Vector<T, N>::dot(const Vector& v) const  
{
  T result=0;
  for (unsigned i=0; i<N; ++i) 
    {
      result += data[i] * v[i];
    }
  return result;
}
    
///  multiplication with scalar
template< typename T, unsigned N> 
inline Vector<T, N> Vector<T, N>::dot(const T a) const  
{
  Vector result;
  for (unsigned i=0; i<N; ++i) 
    {
      result[i] = a* data[i];
    }
  return result;
}
    
/// cross-product-warning
template< typename T, unsigned N> 
inline Vector<T, N> Vector<T, N>::cross(const Vector& v) const  
{
  std::cout << std::endl <<
    "---------------------------" << std::endl <<
    " cross product not defined for N != 3 " << std::endl <<
    "---------------------------" << std::endl;
  assert(false);
  Vector dummy;
  return dummy;
}
        
/// assign addition
template< typename T, unsigned N> 
inline Vector<T, N> & Vector<T, N>::operator += (const Vector& v) 
{
  for (unsigned i=0; i<N; ++i) 
    {
      data[i] += v[i];
    }
  return *this;
}


template< typename T, unsigned N> 
inline Vector<T, N> Vector<T, N>::unit_vec()
{ 
  return *this / mag() ;
}

///////////////////////////////////////////////


// make real vector complex
template< typename T, unsigned N> 
Vector<std::complex<T>, N> Vector<T, N>::make_complex() const
{
  Vector<std::complex<T>, N> result;
  for(unsigned i=0; i<N; ++i)
    result[i] = std::complex<T>(data[i]);

  return result;
}

// make complex vector real
template< typename T, unsigned N> 
Vector<double, N> Vector<T, N>::abs() const
{
  Vector<double, N> result;
  for(unsigned i = 0; i< N; ++i)
    result[i] = std::abs(data[i]);
  
  return result;
}


/// 3D cross product
template<>
inline Vector<double, 3> Vector<double, 3>::cross(const Vector<double, 3>& v) const    
{
    Vector<double, 3> result;
    result[0] = data[1]*v[2]-v[1]*data[2];
    result[1] = data[2]*v[0]-v[2]*data[0];
    result[2] = data[0]*v[1]-v[0]*data[1];
    
    return result;
}



// global methods --> symbols for calculations and printing


template< typename T, unsigned N>  /// Vector * Vector
inline double operator * (const Vector<T, N>& a , const Vector<double, 3>& b ) /// scalar product
{
    return a.dot(b);
}

template< typename T, unsigned N>  /// Vector * scalar
inline Vector<T,N> operator * (const Vector<T, N> & v, const double a)
{
    return v.dot(a);
}

template< typename T, unsigned N>  /// scalar * Vector
inline Vector<T,N> operator * (const double a, const Vector<T, N> & v)
{
    return v.dot(a);
}

template< typename T, unsigned N>  /// Vector / scalar
inline Vector<T,N> operator / (const Vector<T, N> & v, const double a)
{
    return v.dot(1/a);
}

template< typename T, unsigned N>  /// cross product --> Vector % Vector
inline Vector<T,N> operator % (const Vector<T,N> & a, const Vector<T,N> & b)  // cross-product
{
    return a.cross(b);
}


//! \brief output stream used on vector object
/*! @param os output stream 
 @param v vector  */
template< typename T, unsigned N>  /// print Vector
inline std::ostream & operator << (std::ostream & os, const Vector<T,N> & v)
{
    os << "(" ;
    for (unsigned i=0; i<(N-1); ++i) {
        os << v[i] << " , ";
    }
    os << v[N-1] << ")"; 
    return os;
}

