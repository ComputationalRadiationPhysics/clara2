# Copyright 2014-2016 Richard Pausch
#
# This file is part of Clara 2.
#
# Clara 2 is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Clara 2 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Clara 2.
# If not, see <http://www.gnu.org/licenses/>.
#

# make the detector class correctly

CC = g++
CFLAGS = -Wall -O3 -c


all: libDetector.a 

libDetector.a: detector_e_field.o  detector_dft.o detector_fft.o 
	rm -f libDetector.a
	ar cr libDetector.a detector_e_field.o detector_dft.o detector_fft.o 

#copy: libDetector.a
#	cp libDetector.a ..


# E_field detector:
detector_e_field.o: detector_e_field.cpp detector_e_field.hpp vector.hpp utilities.hpp large_index_storage.hpp
	$(CC) $(CFLAGS) detector_e_field.cpp

# Discrete Fourier Transformation detectors:
detector_dft.o: detector_dft.cpp detector_dft.hpp vector.hpp utilities.hpp
	$(CC) $(CFLAGS) detector_dft.cpp

# Fast Fourier Transformation detectors:
detector_fft.o: detector_fft.cpp detector_fft.hpp vector.hpp utilities.hpp fft_ned.hpp
	$(CC) $(CFLAGS) detector_fft.cpp



# How do I include the detector.hpp <-- vector.hpp dependency


clean:
	rm -f *o
	rm -f libDetector.a
