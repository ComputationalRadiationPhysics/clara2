# Copyright 2014-2018 Richard Pausch, Alexander Koehler
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

CXX ?= g++
CFLAGS = -Wall -O3
# compile against fftw library
# http://www.fftw.org/ - an open source FFT library under GPL 2.0 license
CFFT = -lfftw3 -lm
COBJ = -c
OMP = -fopenmp
ZIP = -lz

MPIFLAG= -D__PARALLEL_SETTING__=1
ARRAYFLAG= -D__PARALLEL_SETTING__=2


MPICC=mpic++


CFLAGFORTRAN = #-DCHECK_FOR_FORTRAN_ERROR
# in case one need to check for FORTRAN output errors



all: MPI process

subsystem:
	$(MAKE) -C ./include/


MPI: subsystem MPI_main all_directions.o ./include/libDetector.a single_direction.o
	$(MPICC) $(CFLAGS) $(CFFT) $(OMP) $(ZIP) main.o all_directions.o single_direction.o ./include/libDetector.a ./include/fileExists.o -o executable

ARRAY: subsystem ARRAY_main all_directions.o ./include/libDetector.a single_direction.o
	$(CXX) $(CFLAGS) $(CFFT) $(OMP) $(ZIP) main.o all_directions.o single_direction.o ./include/libDetector.a ./include/fileExists.o -o executable


MPI_main: main.cpp all_directions.hpp
	$(MPICC) $(CFLAGS) $(CFFT) $(MPIFLAG) $(COBJ) $(OMP) $(ZIP) main.cpp

ARRAY_main:  main.cpp all_directions.hpp
	$(CXX) $(CFLAGS) $(CFFT) $(ARRAYFLAG) $(COBJ) $(OMP) $(ZIP) main.cpp

# main routine:
single_direction.o: single_direction.hpp single_direction.cpp ./include/detector_e_field.hpp ./include/detector_dft.hpp \
             ./include/detector_fft.hpp ./include/vector.hpp ./include/import_from_file.hpp ./include/discrete.hpp \
             ./run_through_data.hpp ./include/load_txt.hpp ./include/interpolation.hpp \
             ./include/interpolation.tpp  ./include/fileExists.hpp ./include/fileExists.cpp
	$(CXX) $(CFLAGS) $(CFFT) $(COBJ) $(CFLAGFORTRAN) $(OMP) -I./include/ single_direction.cpp

all_directions.o: all_directions.cpp all_directions.hpp single_direction.hpp ./include/vector.hpp settings.hpp \
			 setFilename.hpp ./include/input_output.hpp
	$(CXX) $(CFLAGS) $(CFFT) $(COBJ) $(OMP) $(ZIP) -I./include/ all_directions.cpp


# How do I include the detector.hpp <-- vector.hpp dependency


process: process_data.cpp settings.hpp setFilename.hpp ./include/input_output.hpp
	$(CXX) $(CFLAGS) $(ZIP) process_data.cpp -o process_data

clean:
	rm -f *o
	rm -f executable
	rm -f process_data
	$(MAKE) -C ./include/ clean
