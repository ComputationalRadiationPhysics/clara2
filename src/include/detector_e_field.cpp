/**
 * Copyright 2014-2018 Richard Pausch
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

#include "physics_units.hpp"
#include "utilities.hpp"


/** constructor for electric field detector (at position in space)
  *
  * @param detector  = R_vec with observation position
  * @param delta_t   = time step width
  * @param N_sig     = maximum number of field entries to the detector
  * @param start_sig = start time at which the electric field should
  *                    be "recorded"
  */
inline Detector_e_field::Detector_e_field(R_vec detector,
                                          double delta_t,
                                          unsigned N_sig,
                                          unsigned start_sig)
  : delta_t(delta_t),
    signal(N_sig, start_sig),
    detector(detector),
    counter(0)
{ }

/** compute integer index in which the retarded electric field (signal)
  * should be written to
  *
  * @param t_signal = time at which the electric field arrives at
  *                   detector position
  * @return unsigned int index in which the time fits best in
  *         the signal class
  */
inline unsigned Detector_e_field::delay_index(double t_signal)
{
  return unsigned(t_signal / delta_t) +1 ;
}


/** compute the time at which the electric field emitted by
  * a charged particle at position r and time t arrives at the
  * detector position (retarded time)
  *
  * @param r = R_vec position of particle
  * @param t = double current time (of the particle)
  * @return signal arrival time (retarded time)
  */
inline double Detector_e_field::t_signal(R_vec r,
                                         double t)
{
  R_vec verbindung = r - detector;  // (? do I use this vector ?)
  double R = verbindung.mag();
  return t + R / phy::c;
}


/** returns the counter of double counts minus no counts
  *
  * @return value of index counter
  */
inline int Detector_e_field::count()
{
  return counter;
}

/* ISSUE #74 - clean up interface */
/** compute electric field value(s) at detector position
  * and store the field values using two time steps
  * (subscript 0 and 1)
  *
  * @param r_0 = particle position at t_part_0
  * @param r_1 = particle position at t_part_0 + delta_t
  * @param p_0 = particle momentum at t_part_0
  * @param p_1 = particle momentum at t_part_0 + delta_t
  * @param dot_p_0 = particle change in momentum at t_part_0
  * @param dot_p_1 = particle change in momentum at t_part_0 + delta_t
  * @param beta_0 = particle beta (v/c) at t_part_0
  * @param beta_1 = particle beta (v/c) at t_part_0 + delta_t
  * @param gamma_0 = particle relativistic gamma at t_part_0
  * @param gamma_1 = particle relativistic gamma at t_part_0 + delta_t
  * @param dot_gamma_0 = particle change in gamma at t_part_0
  * @param dot_gamma_1 = particle change in gamma at t_part_0 + delta_t
  * @param t_part_0 = time at 'first' time step
  */
void Detector_e_field::place(const R_vec r_0,
                             const R_vec r_1,
                             const R_vec p_0,
                             const R_vec p_1,
                             const R_vec dot_p_0,
                             const R_vec dot_p_1,
                             const R_vec beta_0,
                             const R_vec beta_1,
                             const double gamma_0,
                             const double gamma_1,
                             const double dot_gamma_0,
                             const double dot_gamma_1,
                             const double t_part_0)
{
  /* compute electric field (signal) arrival times at detector
   * and map this to a discrete index entry in the detector electric
   * field array */
  double t_signal_0 = t_signal(r_0, t_part_0);
  double t_signal_1 = t_signal(r_1, t_part_0+delta_t);

  unsigned delay_index_0 = delay_index(t_signal_0);
  unsigned delay_index_1 = delay_index(t_signal_1);

  /* ISSUE #78 - time index might be off by one */
  /* in case both time step would be closest to two neighboring
   * array entries interpolate between both time step to fill one
   * array entry */
  if (delay_index_0 == delay_index_1 -1)
  {
    /* time difference between discrete index and retarded time */
    double t_prim =  delta_t * delay_index_0 - t_signal_0;
    /* interpolate between time step 0 and 1 */
    R_vec r_prim = interpol(r_0, r_1, t_prim);
    R_vec p_prim = interpol(p_0, p_1, t_prim);
    R_vec dot_p_prim = interpol(dot_p_0, dot_p_1, t_prim);
    double gamma_prim = interpol(gamma_0, gamma_1, t_prim);
    double dot_gamma_prim = interpol(dot_gamma_0, dot_gamma_1, t_prim);
    R_vec beta_prim = interpol(beta_0, beta_1, t_prim);
    /* calculate beta * gamma at time t_prim (+ t_signal_0) */
    R_vec beta_dot_times_gamma_prim = (dot_p_prim*(1./(phy::m_e*phy::c))
                                       - dot_gamma_prim*beta_prim);
    /* calculate distance between particle and detector position and derive
     * necessary quantities */
    R_vec distance_prim = r_prim - detector;
    R_vec e_R_prim = distance_prim.unit_vec();
    double R_prim = distance_prim.mag();

    /* ISSUE #78 - time index might be off by one */
    /* write electric field into electric field array (signal) */
    signal[delay_index_0] = Lienard_Wiechert(e_R_prim, p_prim,
                                             beta_dot_times_gamma_prim,
                                             gamma_prim, R_prim);
  }

  /* both time steps would be closest to the same time index
   * do not calculate the electric field */
  else if (delay_index_0 == delay_index_1)
  {
    counter--; /* reduce counter by one */
  }

  /* ISSUE #78 - time index might be off by one */
  /* the two time step cover two three discrete time steps
   * do multiple interpolations */
  else if (delay_index_0 +2 == delay_index_1 )
  {
    counter++ ;
    /* time difference between discrete index and retarded time */
    double t_prim1 =  delta_t * delay_index_0 - t_signal_0;

    /* ---- first time index ---- */
    /* interpolate between time step 0 and 1  for first index */
    R_vec r_prim1 = interpol(r_0, r_1, t_prim1);
    R_vec p_prim1 = interpol(p_0, p_1, t_prim1);
    R_vec dot_p_prim1 = interpol(dot_p_0, dot_p_1, t_prim1);
    double gamma_prim1 = interpol(gamma_0, gamma_1, t_prim1);
    double dot_gamma_prim1 = interpol(dot_gamma_0, dot_gamma_1, t_prim1);
    R_vec beta_prim1 = interpol(beta_0, beta_1, t_prim1);
    /* calculate beta * gamma at time t_prim1 (+ t_signal_0) */
    R_vec beta_dot_times_gamma_prim1 = (dot_p_prim1*(1./(phy::m_e*phy::c))
                                        - dot_gamma_prim1*beta_prim1);
    /* calculate distance between particle and detector position and derive
     * necessary quantities */
    R_vec distance_prim1 = r_prim1 - detector;
    R_vec e_R_prim1 = distance_prim1.unit_vec();
    double R_prim1 = distance_prim1.mag();

    /* ISSUE #78 - time index might be off by one */
    /* write electric field into electric field array (signal) */
    signal[delay_index_0] = Lienard_Wiechert(e_R_prim1, p_prim1,
                                             beta_dot_times_gamma_prim1,
                                             gamma_prim1, R_prim1);

     /* ---- second time index ---- */
    /* interpolate between time step 0 and 1  for second index */
    double t_prim2 =  delta_t * (delay_index_0+1) - t_signal_0;
    R_vec r_prim2 = interpol(r_0, r_1, t_prim2);
    R_vec p_prim2 = interpol(p_0, p_1, t_prim2);
    R_vec dot_p_prim2 = interpol(dot_p_0, dot_p_1, t_prim2);
    double gamma_prim2 = interpol(gamma_0, gamma_1, t_prim2);
    double dot_gamma_prim2 = interpol(dot_gamma_0, dot_gamma_1, t_prim2);
    R_vec beta_prim2 = interpol(beta_0, beta_1, t_prim2);
    /* calculate beta * gamma at time t_prim2 (+ t_signal_0) */
    R_vec beta_dot_times_gamma_prim2 = (dot_p_prim2*(1./(phy::m_e*phy::c))
                                        - dot_gamma_prim2*beta_prim2);
    /* calculate distance between particle and detector position and derive
     * necessary quantities */
    R_vec distance_prim2 = r_prim2 - detector;
    R_vec e_R_prim2 = distance_prim2.unit_vec();
    double R_prim2 = distance_prim2.mag();

    /* ISSUE #78 - time index might be off by one */
    /* write electric field into electric field array (signal) */
    signal[delay_index_0+1] = Lienard_Wiechert(e_R_prim2, p_prim2,
                                               beta_dot_times_gamma_prim2,
                                               gamma_prim2, R_prim2);
  }
  /* in case none of the above occurs - throw an error
   * (if more discrete time steps are in between the linear interpolation
   * would hide the electric field dynamic due to a bad temporal resolution) */
  else
  {
    std::cout << "unknown error " << std::endl;
    assert(false);
  }
}


/** simple interpolation between two values (f.e. location, speed)
  * at two time step with time difference delta_t
  *
  * @param x_0   = value at start point
  * @param x_1   = value at endpoint (one time step later)
  * @param t     = time as interpolation parameter (0 < t < delta_t)
  */
template<typename V>
V Detector_e_field::interpol(V x_0,
                             V x_1,
                             double t)
// interpolation between 2 points
{
  return x_0 + ((x_1 - x_0)*(t/delta_t));
}


/** calculates the electric field $\vec E$ based on Lienard Wiechert potential
  * @param e_R      = unit vector pointing from the electron to the detector
  * @param beta     = beta of the electron
  *                   $ \vec \beta = \frac{\vec v}{c} $
  * @param beta_dot = $ \operatorname{\frac{d}{dt}} \vec \beta $
  * @param gamma    = gamma value of the electron
  *                   $ gamma = \sqrt{\frac{1}{1-\vec \beta^2} } $
  * @param R       = distance between electron and detector
  * @return electric field at detector position
  */
R_vec Detector_e_field::Lienard_Wiechert(const R_vec& e_R,
                                         const R_vec& p,
                                         const R_vec beta_dot_times_gamma,
                                         const double& gamma,
                                         const double& R)
{
  R_vec p_o = p*(1/(phy::c*phy::m_e));

  return (1./(4.*M_PI*phy::epsilon_0))
         * ( phy::q * e_R % ((gamma*e_R-p_o) % beta_dot_times_gamma)*gamma
             / (phy::c*R*util::cube(gamma-p_o*e_R)) ) ;
}
