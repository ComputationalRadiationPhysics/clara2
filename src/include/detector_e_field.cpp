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




#include "detector_e_field.hpp"

//! \brief calculates retardated signal at detector
/*! @param r_0      = position one time step before last step
 @param r_1      = position at last time step
 @param p_0      = momentum one time step before last step
 @param p_1      = momentum at last time step
 @param dot_p_0  = dp/dt one time step before last step
 @param dot_p_1  = dp/dt at last time step
 @param beta_0   = v/c = beta one time step befor last step
 @param beta_1   = v/c = beta at last time step
 @param gamma_0  = gamma one time step befor last step
 @param gamma_1  = gamma at last time step
 @param dot_gamma_0 = d gamma/dt one time step befor last step
 @param dot_gamma_1 = d gamma/dt at last time step
 @param t_part_0 = time at the electron one time step before last step
 */ 

inline Detector_e_field::Detector_e_field(R_vec detector, double delta_t, unsigned N_sig, 
				   unsigned start_sig )
  : delta_t(delta_t), signal(N_sig, start_sig), detector(detector),
    counter(0) 
{ }
	
inline unsigned Detector_e_field::delay_index(double t_signal)
{
  return unsigned(t_signal / delta_t) +1 ;
}
    

inline double Detector_e_field::t_signal(R_vec r, double t)
{
  R_vec verbindung = r - detector;  // (? do I use this vector ?)
  double R = verbindung.mag();
  return t + R / phy::c;
}



inline int Detector_e_field::count()
{ 
  return counter; 
}



void Detector_e_field::place(const R_vec r_0, const R_vec r_1, 
						const R_vec p_0, const R_vec p_1, 
						const R_vec dot_p_0, const R_vec dot_p_1,
						const R_vec beta_0, const R_vec beta_1,
						const double gamma_0, const double gamma_1,
						const double dot_gamma_0, const double dot_gamma_1,
						const double t_part_0)
{
    // signal arrival times discrete     
    double t_signal_0 = t_signal(r_0, t_part_0);
    double t_signal_1 = t_signal(r_1, t_part_0+delta_t);
        
    unsigned delay_index_0 = delay_index(t_signal_0);
    unsigned delay_index_1 = delay_index(t_signal_1);
	    
    if (delay_index_0 == delay_index_1 -1)
    {
        // preparation for Lienard Wiechert
        double t_prim =  delta_t * delay_index_0 - t_signal_0;
        R_vec r_prim = interpol(r_0, r_1, t_prim);
        R_vec p_prim = interpol(p_0, p_1, t_prim);
        R_vec dot_p_prim = interpol(dot_p_0, dot_p_1, t_prim);
		
        double gamma_prim = interpol(gamma_0, gamma_1, t_prim);
        double dot_gamma_prim = interpol(dot_gamma_0, dot_gamma_1, t_prim);
        R_vec beta_prim = interpol(beta_0, beta_1, t_prim);
        R_vec beta_dot_times_gamma_prim = (dot_p_prim*(1./(phy::m_e*phy::c))
										   - dot_gamma_prim*beta_prim);
		
        R_vec verbindung_prim = r_prim - detector;
        R_vec e_R_prim = verbindung_prim.unit_vec();
        double R_prim = verbindung_prim.mag();
		
		
        signal[delay_index_0] = Lienard_Wiechert(e_R_prim, p_prim, 
												 beta_dot_times_gamma_prim,
												 gamma_prim, R_prim);        
    }
    
    else if (delay_index_0 == delay_index_1)
    {
        counter-- ;
    }
    
    else if (delay_index_0 +2 == delay_index_1 )
    {
        counter++ ;
        // preparation for Lienard Wiechert
        double t_prim1 =  delta_t * delay_index_0 - t_signal_0;
        R_vec r_prim1 = interpol(r_0, r_1, t_prim1);
        R_vec p_prim1 = interpol(p_0, p_1, t_prim1);
        R_vec dot_p_prim1 = interpol(dot_p_0, dot_p_1, t_prim1);
        
        double gamma_prim1 = interpol(gamma_0, gamma_1, t_prim1);
        double dot_gamma_prim1 = interpol(dot_gamma_0, dot_gamma_1, t_prim1);
        R_vec beta_prim1 = interpol(beta_0, beta_1, t_prim1);
        R_vec beta_dot_times_gamma_prim1 = (dot_p_prim1*(1./(phy::m_e*phy::c))
											- dot_gamma_prim1*beta_prim1);
		
        R_vec verbindung_prim1 = r_prim1 - detector;
        R_vec e_R_prim1 = verbindung_prim1.unit_vec();
        double R_prim1 = verbindung_prim1.mag();
        
        // ---------------------------
		
        double t_prim2 =  delta_t * (delay_index_0+1) - t_signal_0;
        R_vec r_prim2 = interpol(r_0, r_1, t_prim2);
        R_vec p_prim2 = interpol(p_0, p_1, t_prim2);
        R_vec dot_p_prim2 = interpol(dot_p_0, dot_p_1, t_prim2);
		
        double gamma_prim2 = interpol(gamma_0, gamma_1, t_prim2);
        double dot_gamma_prim2 = interpol(dot_gamma_0, dot_gamma_1, t_prim2);
        R_vec beta_prim2 = interpol(beta_0, beta_1, t_prim2);
        R_vec beta_dot_times_gamma_prim2 = (dot_p_prim2*(1./(phy::m_e*phy::c))
											- dot_gamma_prim2*beta_prim2);
		
        R_vec verbindung_prim2 = r_prim2 - detector;
        R_vec e_R_prim2 = verbindung_prim2.unit_vec();
        double R_prim2 = verbindung_prim2.mag();
        
		
        signal[delay_index_0] = Lienard_Wiechert(e_R_prim1, p_prim1, 
												 beta_dot_times_gamma_prim1, 
												 gamma_prim1, R_prim1); 
        signal[delay_index_0+1] = Lienard_Wiechert(e_R_prim2, p_prim2, 
												   beta_dot_times_gamma_prim2, 
												   gamma_prim2, R_prim2);
	}
    
    else
    {
        std::cout << "Unbekannter Fehler " << std::endl;
        assert(false);
    }    
}


//! \brief simple interpolation between two values (f.e. location, speed)
/*! @param r_0   = value at startpoint
 @param r_1   = value at endpoint (one timestep later)
 @param t     = time as interpolation parameter (0 < t < delta_t)
 */ 
template<typename V>
V Detector_e_field::interpol(V r_0, V r_1, double t) 
// interpolation between 2 points
{
    return r_0 + ((r_1 - r_0)*(t/delta_t));
}





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
R_vec Detector_e_field::Lienard_Wiechert(const R_vec& e_R, const R_vec& p, 
					const R_vec beta_dot_times_gamma, const double& gamma, 
					const double& R)
{  
    R_vec p_o = p*(1/(phy::c*phy::m_e));
	
    return (1./(4.*M_PI*phy::epsilon_0))*
	( phy::q * e_R % ((gamma*e_R-p_o) % beta_dot_times_gamma)*gamma / 
	 (phy::c*R*util::cube(gamma-p_o*e_R)) ) ; 
}



