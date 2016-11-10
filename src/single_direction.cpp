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


#include "single_direction.hpp"

#include "interpolation.hpp"
#include "discrete.hpp"
#include "run_through_data.hpp"


/**
 * calculates a single spectra for only one trace and one direction
 *
 * @param data pointer to trajectory data
 * @param linenumber number of data points
 * @param all_omega pointer to frequency values
 * @param all_spectrum pointer to memory for spectra
 * @param N_all_spec maximum number of spectral data allocated
 * @param theta_offset offset of angle theta (used to set angle)
 * @param phi_offset offset of angle phi (used to set angle)
 **/
int single_direction(const one_line* data,
                     const unsigned int linenumber,
                     const double* all_omega,
                     double* all_spectrum,
                     const unsigned N_all_spec,
                     const double theta_offset,
                     const double phi_offset)
{

  /* ----------------------- detectors ------------------------ */
  /* creating the "looking vector" for the "detector class object"
   * that "observes" the emitted radiation in the far field */

  /* setup observation direction: */
  R_vec looking_vector;
  double angle_theta;
  double angle_phi;

  angle_theta = theta_offset/180.0 * M_PI;

  angle_phi =  phi_offset/180.0 * M_PI;
  
  looking_vector = R_vec(std::cos(angle_theta) ,
                         std::sin(angle_theta) * std::cos(angle_phi)  ,
                         std::sin(angle_theta) * std::sin(angle_phi) );


  /* -------- FFT ------------- */
  /* this performs the radiation calculation
   * (for a single trace and a single direction)
   * using the FFT algorithm */


  /* create memory for detectors */
  Detector_fft* detector_fft;

  /* create FFT detectors */
  detector_fft = new Detector_fft(looking_vector, linenumber);



  /* convert trajectory data to meaning full values */
  /* TO DO: FUNCTION CALLED FOR EACH DIRECTION - ISSUE #14 */
  run_through_data(data, linenumber, detector_fft);

  detector_fft->calc_spectrum();

  
  
  /* this switches off the computation using a DFT method - ISSUE #13 */
#if 0

  /* ------------------ DFT  ------------------- */
  /* right now, this is just a debugging method for the FFT
   * it calculates the radiation for the exact same frequencies
   * as the FFT */


  /* create detector objects for radiation */
  Detector_dft* detector_dft;


  /* create DFT detectors (using omega[i] from FFT) to be on the exact same frequency */
  /* number of (use-full) frequencies from FFT */
  unsigned dummy_N =  detector_fft->half_frequency();

  /* get frequencies from FFT */
  double* dummy_omega = detector_fft->frequency;
  detector_dft = new Detector_dft(looking_vector, dummy_N, dummy_omega);

  /* convert trajectory data to meaningful values */
  /* TO DO: FUNCTION CALLED FOR EACH DIRECTION - ISSUE #14 */
  run_through_data(data, linenumber, detector_dft);

  /* calculate spectra using DFT method */
  detector_dft->calc_spectrum();


#endif
  /* end DFT calculation switch - ISSUE #13 */



  /* ------ interpolate spectra onto requested frequencies ---------- */
  /* this is the "integrated interpolation" found in interpolate.tpp */
  interpolation_int(detector_fft, all_omega, all_spectrum, N_all_spec);
  
  
  /* ----- free memory used for detectors ---- */
  
  /* TO DO: DFT detectors not used - ISSUE #13 */
  /* delete detector_dft; */
  delete detector_fft;
  
  return 0;
}
