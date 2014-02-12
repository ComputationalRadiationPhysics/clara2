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





#include "utilities.hpp"
#include "physics_units.hpp"

#ifndef DISCRETE_RPAUSCH
#define DISCRETE_RPAUSCH



/**
 * \brief: storage class to handle 4 values and calculate derivatives
 * usage Descrete<datatype> 
 */
 
template<typename T>
class Discrete {
    //  friend class Retardation;
public:
    
    //! \brief constructor filling all values
    /*! @param old2 value at time = t-3
        @param old  value at time = t-2
        @param now  value at time = t-1
        @param future value at time = t
        @param h pointer to discrete of time */
    Discrete(T old2, T old, T now, T future, const Discrete<double>* h = 0)
        : old2(old2), old(old), now(now), future(future), h(h) {}
    
    //! \brief constructor without filling data
    /*! @param h difference between time steps */
    Discrete(const Discrete<double>* h) 
        : h(h) {}
    
    //! \brief copy constructor
    Discrete<T>& operator=(const Discrete<T>& copy)
    {
        assert(copy.h == h);
        old2   = copy.old2;
        old    = copy.old;
        now    = copy.now;
        future = copy.future;
        
        return *this;
    } 
    
    //! \brief setting new future value and moving data down (now --> old, ...)
    //! @param next new future value 
    void next(T next)
    {
        old2 = old;
        old = now;
        now = future;
        future = next;
    }
    
    //! \brief returning derivative for time = t - 2
    T dot_old() const
    {return (now - old2)/(h->get_now() - h->get_old2()); }
    
    //! \brief returning derivative for time = t - 1
    T dot_now() const
    {return (future - old)/(h->get_future() - h->get_old()); }
    
    //! \brief return value for time = t - 3 (old2)
    T get_old2() const
    {return old2; }
    
    //! \brief return value for time = t - 2 (old)
    T get_old() const 
    {return old; }
    
    //! \brief return value for time = t - 1 (now)
    T get_now() const
    {return now; }
    
    //! \brief return value for time = t - 0 (future)
    T get_future() const
    {return future; }

    //! \brief return value for time = t - 0 (future)
    T get_delta_old() const
    {return (get_now()-get_old()); }
    
    
private:
    T old2, old, now, future;
    const Discrete<double>* h;  
};


/*!
 * \brief a class provinding additional methods for relativistic physics\n
 * returns Discrete objects for gamma and beta\n
 * and also provides functors for calculating single values of gamma and beta\n
 * Energie = sqrt(p^2*c^2+m_0^2*c^4) = GAMMA*m
 * \vec beta = (\vec v)/c  with \vec v = SPEED = (\vec p)/m(v) 
 * = (\vec p)/(m_0 * gamma)  
 */
class More_discrete
{
public:
    //! \brief constructor 
    /*! @param Det is a reference to a detector class from which one gets:
        -> stepwidth: the length of the timestep between to steps */
    More_discrete(Discrete<double>* h)
        : stepwidth(h) {}
    
    
    //! \brief returns a Discrete object for gamma
    /*! @param p is a Discrete object of the momentum */ 
    Discrete<double> momentum_to_gamma(Discrete<R_vec>& p)
    {
        return Discrete<double>(gamma(p.get_old2()), 
                                gamma(p.get_old()), 
                                gamma(p.get_now()), 
                                gamma(p.get_future()), 
                                stepwidth);
    }
    
    //! \brief returns a Discrete object for beta
    /*! @param p is a Discrete object of the momentum 
        @param gamma is a Discrete object of gamma --> better to calculate here?
     */     
    Discrete<R_vec> momentum_to_beta(Discrete<R_vec> p, Discrete<double> gamma)
    {
        return Discrete<R_vec>(beta(p.get_old2(), gamma.get_old2()),
                               beta(p.get_old(), gamma.get_old()),
                               beta(p.get_now(), gamma.get_now()), 
                               beta(p.get_future(), gamma.get_future()),
                               stepwidth);    
    } 
    
    //! \brief functor to calculate a single gamma value
    //! @param p momentum
    double gamma(R_vec p){
        return sqrt(util::square<R_vec, double>(p*phy::c)+
                    util::square(phy::m_e*util::square(phy::c)))/(phy::m_e*
                    util::square(phy::c));
    }

    //! \brief functor to calculate a single beta value
    /*! @param p momentum
        @param gamma gamma */
    R_vec beta(R_vec p, double gamma)
    {
        return p*(1.0/(phy::c*phy::m_e*gamma));
    }
    
private:
    const Discrete<double>* stepwidth;
    
};


#endif

