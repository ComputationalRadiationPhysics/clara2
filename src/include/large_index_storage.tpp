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



#include "large_index_storage.hpp"



template< typename T>
inline Large_index_storage<T>::Large_index_storage(const unsigned N, const unsigned start)
  : N(N), start(start), signal(new T[N]) 
{ }
    
template<typename T>
inline Large_index_storage<T>:: ~Large_index_storage()
{
  delete[] signal;
}

//! \brief direct access to data, returns reference 
/*! \param i is equivalent to x[i] */
template<typename T>
inline T& Large_index_storage<T>::operator[](unsigned i)
{
  if (!((i>=start) && (i<start+N))) {
    std::cout << "Index falsch: " << i << " -> " 
	      << start << " - " << start+N<<std::endl;
    assert(false); // aborts a wrong access to data
  }
  return signal[i-start]; // returns the data 
}

//! \brief returns stored data and zero if not in [start, start + N -1]  
/*! \param i is equivalent to x[i] */
template<typename T>
inline T Large_index_storage<T>::operator()(unsigned i)
{
  if((i>= start) && (i< start+N))
    return signal[i-start]; // returns data
  else 
    return T(0);            // returns zero (or equivalent)
}


