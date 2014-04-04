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



#ifndef INTERPOLATION_RPAUSCH
#define INTERPOLATION_RPAUSCH

#include "detector_fft.hpp"

template <typename X, typename Y>
void interpolation(const X* x_old, const Y* y_old, const unsigned N_old, 
		   const X* x_new,       Y* y_new, const unsigned N_new);



template <typename X, typename Y>
void interpolation_on(const X* x_old, const Y* y_old, const unsigned N_old, 
		      const X* x_new,       Y* y_new, const unsigned N_new);

template <typename X, typename Y>
void interpolation_on(const Detector_fft* fft, 
		      const X* x_new,       Y* y_new, const unsigned N_new);


void interpolation_int(Detector_fft* fft, 
		       const double* x_new,       double* y_new, const unsigned N_new);



#include "interpolation.tpp"

#endif
