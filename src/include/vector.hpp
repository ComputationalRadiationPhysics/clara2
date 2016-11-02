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
#include <cmath>
#include <cassert>
#include <complex>

#pragma once

/**
 * \brief  Vector class by Richard
 *
 *  This class provides operations like: +, -, *, %, /, =, +=, [], mag() \n
 *  usage: Vector<datatype, size> \n
 *  short for Vector<double, 3> --> R_vec
 */

template< typename T, unsigned N>
class Vector{

public:
  // constructor, destructor:
  /// basic constructor setting everything to zero
  inline Vector();

  //! \brief special constructor for 3D
  /*! @param x first coordinate
      @param y second coordinate
      @param z third coordinate */
  inline Vector(double x,
                double y=0.,
                double z=0.);

  /// copy constructor --> not necessary

  /// destructor
  inline ~Vector() {}  // not necessary


  // Getters and Setters:
  /// access data
  inline T & operator[](unsigned i);

  /// access data
  const inline T & operator[](unsigned i) const;



  // Calculations:

  /// addition
  inline Vector operator+ (const  Vector& v) const;

  /// subtraction
  inline Vector operator- (const Vector& v) const;

  /// magnitude
  inline T mag() const;

  /// scalar product
  inline T dot(const Vector& v) const;

  ///  multiplication with scalar
  inline Vector dot(const T a) const;

  /// cross-product-warning
  inline Vector cross(const Vector& v) const;

  /// assign addition
  inline Vector & operator += (const Vector& v);

  inline Vector unit_vec();


  ////////////

  /// make real vector complex
  Vector<std::complex<T>, N> make_complex() const;

  /// make complex vector real (absolute value of complex number)
  Vector<double, N> abs() const;



private:
  /// stored data
  T data[N];

};


// Member function specialization:


/// 3D cross product
template<>
inline Vector<double, 3> Vector<double, 3>::cross(const Vector<double, 3>& v) const;



// global methods --> symbols for calculations and printing

template< typename T, unsigned N>  /// Vector * Vector --> scalar product
inline double operator * (const Vector<T, N>& a ,
                          const Vector<double, 3>& b );

template< typename T, unsigned N>  /// Vector * scalar
inline Vector<T,N> operator * (const Vector<T, N> & v,
                               const double a);

template< typename T, unsigned N>  /// scalar * Vector
inline Vector<T,N> operator * (const double a,
                               const Vector<T, N> & v);

template< typename T, unsigned N>  /// Vector / scalar
inline Vector<T,N> operator / (const Vector<T, N> & v,
                               const double a);

template< typename T, unsigned N>  /// cross product --> Vector % Vector
inline Vector<T,N> operator % (const Vector<T,N> & a,
                               const Vector<T,N> & b);


//! \brief output stream used on vector object
/*! @param os output stream
    @param v vector  */
template< typename T, unsigned N>  /// print Vector
inline std::ostream & operator << (std::ostream & os,
                                   const Vector<T,N> & v);


/*! \var typedef Vector<double, 3> R_vec
  \brief A type definition for a 3D double vector.

  Because in physics this is the most widely used vector, there is a special
  typedef for it.
*/

#include "vector.tpp"

typedef Vector<double, 3> R_vec;
