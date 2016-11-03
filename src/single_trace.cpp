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


#include "single_trace.hpp"

#include "interpolation.hpp"
#include "discrete.hpp"
#include "run_through_data.hpp"


/**
 * calculates a single spectra for only one trace and one direction
 *
 * @param one_line pointer to trajectory data
 * @param linenumber number of data points
 * @param all_omega pointer to frequency values
 * @param all_spectrum pointer to memory for spectra
 * @param N_all_spec maximum number of spectral data allocated
 * @param theta_offset offset of angle theta (used to set angle)
 * @param phi_offset offset of angle phi (used to set angle)
 **/
int single_trace(const one_line* data,
                 const unsigned int linenumber,
                 const double* all_omega,
                 double* all_spectrum,
                 const unsigned N_all_spec,
                 const double theta_offset,
                 const double phi_offset)
{
  ////////////////////////////////////////////////////////////
  /* CHANGE THIS: */
  /* this are the remaining of the old program
   * structure, where several directions were calculated
   * right here - ISSUE #12 */

  const unsigned N_angle_theta = 1;
  const unsigned N_angle_phi   = 1;


  //////////////////////////////////////////////////////////////


  /* ----------------------- detectors ------------------------ */
  /* creating the "looking vector" for the "detector class object"
   * that "observes" the emitted radiation in the far field */

  /* creating different angles: */
  /* TO DO: SINCE ONLY A SINGLE DIRECTION IS CALCULATED IN THIS FUNCTION
   * THE ARRAY SETUP IS OBSOLETE  - ISSUE #12 */
  R_vec looking_vector[N_angle_theta * N_angle_phi];
  double angle_theta[N_angle_theta];
  double angle_phi[N_angle_phi];

  for(unsigned i=0; i< N_angle_theta; ++i)
  {
    angle_theta[i] = (double(i)/double(N_angle_theta)) * M_PI * 0.5 + theta_offset/180.0 * M_PI;
  }

  for(unsigned i=0; i< N_angle_phi; ++i)
  {
    angle_phi[i] = (double(i)/double(N_angle_phi))*1.0 * M_PI + phi_offset/180.0 * M_PI;
  }

  for (unsigned a=0; a<N_angle_theta; ++a)
  {
    for (unsigned b=0; b<N_angle_phi; ++b)
    {
      looking_vector[a*N_angle_phi+b] = R_vec(std::cos(angle_theta[a]) ,
                                              std::sin(angle_theta[a]) * std::cos(angle_phi[b])  ,
                                              std::sin(angle_theta[a]) * std::sin(angle_phi[b]) );
    }
  }


  /* -------- FFT ------------- */
  /* this performs the radiation calculation
   * (for a single trace and a single direction)
   * using the FFT algorithm */


  /* create memory for detectors */
  /* TO DO: SEVERAL DIRECTIONS - ISSUE #12 */
  Detector_fft* detector_fft[N_angle_theta*N_angle_phi];

  /* create FFT detectors */
  for(unsigned i=0; i<N_angle_theta*N_angle_phi; ++i)
  {
    detector_fft[i] = new Detector_fft((looking_vector[i]), linenumber );
  }


  /* convert trajectory data to meaning full values */
  /* TO DO: FUNCTION CALLED FOR EACH DIRECTION - ISSUE #14 */
  run_through_data(data, linenumber, N_angle_theta*N_angle_phi, detector_fft);

  /* calculate spectra FFT */
  for(unsigned k=0; k<N_angle_theta*N_angle_phi; ++k)
  {
    (*detector_fft[k]).calc_spectrum();
  }


  /* this switches off the computation using a DFT method - ISSUE #13 */
#if 0

  /* ------------------ DFT  ------------------- */
  /* right now, this is just a debugging method for the FFT
   * it calculates the radiation for the exact same frequencies
   * as the FFT */


  /* create detector objects for radiation */
  /* TO DO: SEVERAL DIRECTIONS - ISSUE #12 */
  Detector_dft* detector_dft[N_angle_theta * N_angle_phi];


  /* create DFT detectors (using omega[i] from FFT) to be on the exact same frequency */
  for(unsigned i=0; i<N_angle_theta*N_angle_phi ; ++i)
  {
    /* number of (use-full) frequencies from FFT */
    unsigned dummy_N =  (*detector_fft[i]).half_frequency();

    /* get frequencies from FFT */
    double* dummy_omega = (*detector_fft[i]).frequency;
    detector_dft[i] = new Detector_dft(looking_vector[i], dummy_N, dummy_omega);
  }

  /* convert trajectory data to meaningful values */
  /* TO DO: FUNCTION CALLED FOR EACH DIRECTION - ISSUE #14 */
  run_through_data(data, linenumber, N_angle_theta*N_angle_phi, detector_dft);

  /* calculate spectra using DFT method */
  for(unsigned k=0; k< N_angle_theta*N_angle_phi ; ++k)
  {
    (*detector_dft[k]).calc_spectrum();
  }

#endif
  /* end DFT calculation switch - ISSUE #13 */



  /* ------ interpolate spectra onto requested frequencies ---------- */
  /* this is the "integrated interpolation" found in interpolate.tpp */
  interpolation_int(detector_fft[0], all_omega, all_spectrum, N_all_spec);


  /* ----- free memory used for detectors ---- */

  for(unsigned i=0; i<N_angle_theta*N_angle_phi; ++i)
  {
    /* TO DO: DFT detectors not used - ISSUE #13 */
    /* delete detector_dft[i]; */
    delete detector_fft[i];
  }

  return 0;
}
