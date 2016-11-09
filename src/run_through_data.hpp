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
#include "settings.hpp"


inline R_vec beta_times_gamma(R_vec beta)
{
  return std::sqrt(1./(1.-util::square<R_vec, double>(beta))) * beta;
}


template<typename DET>
void run_through_data(const one_line* data,
                      const unsigned linenumber,
                      DET detector)
{
  /* ---------- storing data : comparable to real data stream (not like a file here) --- */

  //time
  Discrete<double> time_fill;
  const Discrete<double> *time_fill_ref = &time_fill;
  //position
  Discrete<R_vec> location( time_fill_ref );
  //momentum
  Discrete<R_vec> momentum( time_fill_ref );

  // decrived quantities
  More_discrete auto_fill( time_fill_ref );
  // beta
  Discrete<R_vec> beta( time_fill_ref );
  // gamma
  Discrete<double> gamma( time_fill_ref );

  /* -------- streaming the data and sending it to the detectors: ---------- */

  namespace in = param::input;
  for(unsigned i=0; i<linenumber; ++i)
  {
    // fill Discrete class with values:
    // time: to seconds
    time_fill.next(double(data[i].intern_data[in::index_time]) * in::convert_time );
    // position: to meter
    location.next( R_vec(data[i].intern_data[in::index_pos_x] * in::convert_pos ,
                         data[i].intern_data[in::index_pos_y] * in::convert_pos ,
                         data[i].intern_data[in::index_pos_z] * in::convert_pos ));
    // momentum: to beta*gamma  --> phy::m_e*beta_times_gamma(beta)*phy:c --> kg*m/s
    const double btom = phy::m_e*phy::c;
    momentum.next( btom*beta_times_gamma(R_vec(data[i].intern_data[in::index_beta_x] * in::convert_beta ,
                                               data[i].intern_data[in::index_beta_y] * in::convert_beta ,
                                               data[i].intern_data[in::index_beta_z] * in::convert_beta ) ));

    gamma.next(auto_fill.gamma(momentum.get_future()));
    beta.next(auto_fill.beta(momentum.get_future(), gamma.get_future()));

    /* if all 4 entries of the Discrete class are filed with values
       the radiation amplitude can be calculated */
    if(i >=3 )
    {
      detector->add_to_spectrum(location.get_old(),
                                beta.get_old(),
                                beta.dot_old(),
                                time_fill.get_old(),
                                time_fill.get_delta_old());
    }

  }

}
