#!/bin/bash

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



#---------------------------------------------------------------
# INFO:
# Running this shell script generates a submit file either
# for a MPI parallelization or a PBS-Array-Job parallelization.
#---------------------------------------------------------------


#------------
# USER SETUP:
#------------

# name of cluster job:
JOBNAME=MPI_lib1

# name of submit file:
SUBMITFILE=submit.sh

# name of queue to submit job to:
QUEUE=laser

# name of file with environment setup:
MODULES2LOAD=clara2_hypnos.modules

# maximum simulation time
WALLTIME=02:00:00

# Job parallel structure:
# for MPI (number of node with (each) number of cores):
NUMBERNODES=4
NUMBERCORES=32
# for Array Jobs:
# set number of parallel tasks to the same as with MPI
ARRAYMAX=$(echo $NUMBERNODES " * " $NUMBERCORES " - 1 " | bc )

#---------------
# END USER SETUP
#---------------


#-------------------------
# RUN TIME USE INTERACTION
#-------------------------

# user interaction required to select parallelization
# present options:
echo "for MPI enter 1"
echo "for PBS-Array jobs enter 2"

# read chosen parallelization method
read DECISION
case  $DECISION in
    1) echo "MPI choosen";;
    2) echo "PBS-Array jobs choosen";;
    *) echo "none choosen"
       exit 1 ;;
esac

#--------------------
# END USE INTERACTION
#--------------------


#-----------------------
# GENERATE SUBMIT SCRIPT
#-----------------------

# generate submit file (here the same for MPI and PBS-array jobs)
echo "#PBS -q "$QUEUE > $SUBMITFILE
echo "#PBS -l walltime="$WALLTIME >> $SUBMITFILE
echo "#PBS -N "$JOBNAME >> $SUBMITFILE
echo "#PBS -m n -M x.lastname@hzdr.de" >> $SUBMITFILE
echo "#PBS -d ." >> $SUBMITFILE
echo "#PBS -e ./error.txt" >> $SUBMITFILE
echo "#PBS -o ./output.txt" >> $SUBMITFILE
echo  "" >> $SUBMITFILE

# in case of a MPI parallelization:
if [ $DECISION -eq 1 ]
then
    echo "#PBS -l nodes="$NUMBERNODES":ppn="$NUMBERCORES >> $SUBMITFILE
else
    # in case of a PBS-array parallelisation
    echo "#PBS -l nodes=1:ppn=1" >> $SUBMITFILE
    echo "#PBS -t 0-"$ARRAYMAX >> $SUBMITFILE
    echo "ARRAYMAX="$ARRAYMAX >> $SUBMITFILE
    echo "export ARRAYMAX" >> $SUBMITFILE
fi

# add newline to submit file and load modules
echo  "" >> $SUBMITFILE
echo  "source " $MODULES2LOAD >> $SUBMITFILE
echo  "" >> $SUBMITFILE

# in case of MPI parallelisation
if [ $DECISION -eq 1 ]
then
    echo -n "mpiexec --prefix $MPIHOME -tag-output --display-map -x LIBRARY_PATH -x LD_LIBRARY_PATH -n " >> $SUBMITFILE
    echo  $(echo $NUMBERNODES " * " $NUMBERCORES | bc )" ./executable  |  tee output.log" >> $SUBMITFILE
else
    # in case of PBS-array parallelization
    echo  " ./executable" >> $SUBMITFILE
fi

# add newline to submit file
echo  "" >> $SUBMITFILE

#---------------------------
# END GENERATE SUBMIT SCRIPT
#---------------------------


#-------------
# COMPILE CODE
#-------------

# load bash environment for compilation
source $MODULES2LOAD

# ISSUE #65 why is make required twice?
make
if [ $DECISION -eq 1 ]
then
    make MPI
else
    make ARRAY
fi

#-----------------
# END COMPILE CODE
#-----------------
