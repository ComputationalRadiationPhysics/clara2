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



#ifndef UTILITIES_RPAUSCH
#define UTILITIES_RPAUSCH

namespace util {

    //goal: to increase readability of code
    
    template<typename A> /// a generic square function
    inline A square(A a )
    {
        return a*a;
    }
    
    template<typename A, typename R> /// a more generic square function
    inline R square(A a )
    {
        return a*a;
    }
    

    template<typename A> /// a generic cube function
    inline A cube(A a)
    {
        return a*a*a;
    }

    template<typename A, typename R> /// a more generic cube function
    inline R cube(A a)
    {
        return a*a*a;
    }
    
}

#endif



