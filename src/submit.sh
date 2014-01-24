# Copyright 2014 Richard Pausch
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

#PBS -q laser
#PBS -l walltime=02:00:00
#PBS -N MPI_lib1
#PBS -m n -M r.pausch@hzdr.de
#PBS -d .
#PBS -e ./error.txt
#PBS -o ./output.txt

#PBS -l nodes=4:ppn=32

mpiexec -n 128 ./executable

