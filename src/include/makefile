# Copyright 2014-2016 Richard Pausch, Alexander Koehler
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
# compile against fftw library
# http://www.fftw.org/ - an open source FFT library under GPL 2.0 license
CFFT = -lfftw3 -lm

all: libDetector.a fileExists.o

libDetector.a: detector_e_field.o  detector_dft.o detector_fft.o
	rm -f libDetector.a
	ar cr libDetector.a detector_e_field.o detector_dft.o detector_fft.o


# E_field detector:
detector_e_field.o: detector_e_field.cpp detector_e_field.hpp vector.hpp utilities.hpp large_index_storage.hpp
	$(CC) $(CFLAGS) detector_e_field.cpp

# Discrete Fourier Transformation detectors:
detector_dft.o: detector_dft.cpp detector_dft.hpp vector.hpp utilities.hpp
	$(CC) $(CFLAGS) detector_dft.cpp

# Fast Fourier Transformation detectors:
detector_fft.o: detector_fft.cpp detector_fft.hpp vector.hpp utilities.hpp ned_fft.hpp ../settings.hpp
	$(CC) $(CFLAGS) $(CFFT) detector_fft.cpp

fileExists.o: fileExists.hpp fileExists.cpp
	$(CC) $(CFLAGS) fileExists.cpp

# How do I include the detector.hpp <-- vector.hpp dependency


clean:
	rm -f *o
	rm -f libDetector.a
