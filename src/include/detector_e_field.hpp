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
#include <cassert>
#include <complex>
#include <cmath>

#include "vector.hpp"
#include "large_index_storage.hpp"
#include "physics_units.hpp"
#include "utilities.hpp"



#ifndef DETECTOR_E_FIELD_RPAUSCH
#define DETECTOR_E_FIELD_RPAUSCH


//! \brief class for a point-like detector storing the signal externaly
class Detector_e_field
{
public:
    //! \brief constructor for a point-like detector
    /*! @param detector = location of the detector 
     @param delta_t  = timestep of odint
	 @param N_sig     = number of datapoints of signal to store
	 @param start_sig = start index of signal
     signal    = pointer to signal struct array (signal at detector)    
     */
    inline Detector_e_field(R_vec detector, double delta_t, unsigned N_sig, 
		     unsigned start_sig );
	

    //! \brief Calculate delay_index from signal time
    //! @param t_signal Time at which the signal arrives at the detector
    inline unsigned delay_index(double t_signal);

    
    //! \brief Calculate time of signal at the detector
    /*! @param r position of moving charge (electron)
	 @param t time at which the particle moved */
    inline double t_signal(R_vec r, double t);

	
    // for details see .cpp file:
    void place(const R_vec r_0, const R_vec r_1, 
               const R_vec p_0, const R_vec p_1, 
               const R_vec dot_p_0, const R_vec dot_p_1,
               const R_vec beta_0, const R_vec beta_1,
               const double gamma_0, const double gamma_1,
               const double dot_gamma_0, const double dot_gamma_1,
               const double t_part_0);
    
    
    /// returns the counter of double counts minus no counts
  inline int count();

	
// data:
    const double delta_t;  /// length of timesteps
    Large_index_storage<R_vec> signal;         /// E_field at detector
	
	
private:
	
    //data
    R_vec detector;        /// location of the detector
    int counter;

	
	//methodes
    
    //! \brief simple interpolation between two values (f.e. location, speed)
    /*! @param r_0   = value at startpoint
     @param r_1   = value at endpoint (one timestep later)
     @param t     = time as interpolation parameter (0 < t < delta_t)
     */ 
    template<typename V>
    V interpol(V r_0, V r_1, double t); // interpolation between 
	// 2 points
    
    /*! \brief calculates Lienard Wiechert Potation (more precise the 
     $\vec E$-field */ 
    /*! @param e_R      = unit vector in the from the electron to the detector
	 @param beta     = beta vector of the electron 
	 $ \vec \beta = \frac{\vec v}{c} $
	 @param beta_dot = $ \operatorname{\frac{d}{dt}} \vec \beta $
	 @param gamma    = gamma value of the electron 
	 $ gamma = \sqrt{\frac{1}{1-\vec \beta^2} } $
	 @ param R       = distance between electron and detector     
     */    
    R_vec Lienard_Wiechert(const R_vec& e_R, const R_vec& p, 
                           const R_vec beta_dot_times_gamma, 
                           const double& gamma, 
                           const double& R);
    
	
	
};


#endif

