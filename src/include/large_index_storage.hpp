/**
 * Copyright 2014 Richard Pausch
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
#include<cassert>
//#include "detector.hpp"

#ifndef LARGE_INDEX_STORAGE_RPAUSCH
#define LARGE_INDEX_STORAGE_RPAUSCH

/**
 * \brief  storage class by Richard 
 *  
 * this class provides a storage systems for large arrays for which
 * there only a few non zero values, which are additionaly close together \n
 * usage: Large_index_storage<datatype> n
 */


template<typename T>
class Large_index_storage
{
public:
// constructor, destructor:
    //! \brief constructor 
    /*! @param N number of posible non zero values (just single segment)
        @param start starting point of array f.e. x[start - x] x >0 --> is not
                defined */
  inline Large_index_storage(const unsigned N, const unsigned start);
    
  inline ~Large_index_storage();

    //! \brief direct access to data, returns reference 
    /*! \param i is equivalent to x[i] */
  inline T& operator[](unsigned i);

    //! \brief returns stored data and zero if not in [start, start + N -1]  
    /*! \param i is equivalent to x[i] */
  inline T operator()(unsigned i);

private:
    const unsigned N, start;
    T* signal;


};

#include "large_index_storage.tpp"


#endif

