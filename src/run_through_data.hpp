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
#include <cstdlib>

#pragma once

#include "discrete.hpp"
#include "import_from_file.hpp"
#include "physics_units.hpp"
#include "vector.hpp"


inline R_vec beta_times_gamma(R_vec beta)
{
  return std::sqrt(1./(1.-util::square<R_vec, double>(beta))) * beta;
}


template<typename DET>
void run_through_data(const one_line* data,
                      const unsigned linenumber,
                      const unsigned N_angle,
                      DET detector)
{
  /* ---------- storing data : comparable to real data stream (not like a file here) --- */

  // USING: SI-UNITS !!!

  //time s --> s
  Discrete<double> time_fill(data[0].intern_data[6] /* *1.E-15 */ ,
                             data[1].intern_data[6] /* *1.E-15 */ ,
                             data[2].intern_data[6] /* *1.E-15 */ ,
                             data[3].intern_data[6] /* *1.E-15 */ );

  //position: in m  --> m
  Discrete<R_vec> location( R_vec(data[0].intern_data[0] /* *1.E-2 */ ,
                                  data[0].intern_data[1] /* *1.E-2 */ ,
                                  data[0].intern_data[2] /* *1.E-2 */ ),
                            R_vec(data[1].intern_data[0] /* *1.E-2 */ ,
                                  data[1].intern_data[1] /* *1.E-2 */ ,
                                  data[1].intern_data[2] /* *1.E-2 */ ),
                            R_vec(data[2].intern_data[0] /* *1.E-2 */ ,
                                  data[2].intern_data[1] /* *1.E-2 */ ,
                                  data[2].intern_data[2] /* *1.E-2 */  ),
                            R_vec(data[3].intern_data[0] /* *1.E-2 */ ,
                                  data[3].intern_data[1] /* *1.E-2 */ ,
                                  data[3].intern_data[2] /* *1.E-2 */  ),
                            &time_fill);

  //momentum: beta*gamma  --> phy::m_e*beta_times_gamma(beta)*phy:c --> kg*m/s
  double btom = phy::m_e*phy::c;
  Discrete<R_vec> momentum( btom*beta_times_gamma(R_vec(data[0].intern_data[3],
                                                        data[0].intern_data[4],
                                                        data[0].intern_data[5]) ),
                            btom*beta_times_gamma(R_vec(data[1].intern_data[3],
                                                        data[1].intern_data[4],
                                                        data[1].intern_data[5]) ),
                            btom*beta_times_gamma(R_vec(data[2].intern_data[3],
                                                        data[2].intern_data[4],
                                                        data[2].intern_data[5]) ),
                            btom*beta_times_gamma(R_vec(data[3].intern_data[3],
                                                        data[3].intern_data[4],
                                                        data[3].intern_data[5]) ),
                            &time_fill);

  More_discrete auto_fill(&time_fill);


  Discrete<R_vec> beta(&time_fill);
  Discrete<double> gamma(&time_fill);

  gamma = auto_fill.momentum_to_gamma(momentum);
  beta  = auto_fill.momentum_to_beta(momentum, gamma);


  /* -------- streaming the data and sending it to the detectors: ---------- */

  for(unsigned i=4; i<linenumber; ++i)
  {
    for(unsigned k=0; k<N_angle; ++k)
    {
      (*detector[k]).add_to_spectrum(location.get_old(),
                                     beta.get_old(),
                                     beta.dot_old(),
                                     time_fill.get_old(),
                                     time_fill.get_delta_old());
    }


    // set new to old:
    time_fill.next(double(data[i].intern_data[6]) /* *1.E-15 */ );

    location.next( R_vec(data[i].intern_data[0] /* *1.E-2 */ ,
                         data[i].intern_data[1] /* *1.E-2 */ ,
                         data[i].intern_data[2] /* *1.E-2 */  ));

    momentum.next( btom*beta_times_gamma(R_vec(data[i].intern_data[3],
                                               data[i].intern_data[4],
                                               data[i].intern_data[5] ) ));

    gamma.next(auto_fill.gamma(momentum.get_future()));

    beta.next(auto_fill.beta(momentum.get_future(), gamma.get_future()));

  }

}
