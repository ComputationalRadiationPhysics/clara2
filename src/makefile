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
CFLAGS = -Wall -O3
COBJ = -c
OMP = -fopenmp
ZIP = -lz

MPIFLAG= -D__PARALLEL_SETTING__=1
ARRAYFLAG= -D__PARALLEL_SETTING__=2


MPICC=mpic++


CFLAGFORTRAN = #-DCHECK_FOR_FORTRAN_ERROR 
# in case one need to check for FORTRAN output errors



all: subsystem MPI process convert_to_matrix

subsystem: 
	$(MAKE) -C ./include/


MPI: MPI_main all_directions.o ./include/libDetector.a single_trace.o
	$(MPICC) $(CFLAGS) $(OMP) $(ZIP) main.o all_directions.o single_trace.o ./include/libDetector.a -o executable

ARRAY: ARRAY_main all_directions.o ./include/libDetector.a single_trace.o
	$(CC) $(CFLAGS) $(OMP) $(ZIP) main.o all_directions.o single_trace.o ./include/libDetector.a -o executable


MPI_main: main.cpp all_directions.hpp
	$(MPICC) $(CFLAGS) $(MPIFLAG) $(COBJ) $(OMP) $(ZIP) main.cpp

ARRAY_main:  main.cpp all_directions.hpp
	$(CC) $(CFLAGS) $(ARRAYFLAG) $(COBJ) $(OMP) $(ZIP) main.cpp

# main routine:
single_trace.o: single_trace.hpp single_trace.cpp ./include/detector_e_field.hpp ./include/detector_dft.hpp \
             ./include/detector_fft.hpp ./include/vector.hpp ./include/import_from_file.hpp ./include/discrete.hpp \
             ./run_through_data.hpp ./include/load_txt.hpp ./include/analytical_solution.hpp ./include/interpolation.hpp \
             ./include/interpolation.tpp
	$(CC) $(CFLAGS) $(COBJ) $(CFLAGFORTRAN) $(OMP) -I./include/ single_trace.cpp

all_directions.o: all_directions.cpp all_directions.hpp single_trace.hpp ./include/vector.hpp ./include/string_manipulation.hpp
	$(CC) $(CFLAGS) $(COBJ) $(OMP) $(ZIP) -I./include/ all_directions.cpp


# How do I include the detector.hpp <-- vector.hpp dependency


process: process_data.cpp
	$(CC) $(CFLAGS) $(ZIP) process_data.cpp -o process_data

convert_to_matrix: convert_to_matrix.cpp
	$(CC) $(CFLAGS) convert_to_matrix.cpp -o convert_to_matrix

clean:
	rm -f *o
	rm executable
	rm process_data
	rm convert_to_matrix
