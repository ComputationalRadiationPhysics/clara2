#!/bin/bash

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

JOBNAME=MPI_lib1
SUBMITFILE=submit.sh
QUEUE=laser

WALLTIME=02:00:00

#Jobstructure:
#for MPI:
NUMBERNODES=4
NUMBERCORES=32
#for Array Jobs:
ARRAYMAX=$(echo $NUMBERNODES " * " $NUMBERCORES " - 1 " | bc )

echo "for MPI enter 1"
echo "for Array jobs enter 2"

read DECISION
case  $DECISION in
    1) echo "MPI choosen";;
    2) echo "Array jobs choosen";;
    *) echo "none choosen"
       exit 1 ;;
esac


echo "#PBS -q "$QUEUE > $SUBMITFILE
echo "#PBS -l walltime="$WALLTIME >> $SUBMITFILE
echo "#PBS -N "$JOBNAME >> $SUBMITFILE
echo "#PBS -m n -M r.pausch@hzdr.de" >> $SUBMITFILE
echo "#PBS -d ." >> $SUBMITFILE
echo "#PBS -e ./error.txt" >> $SUBMITFILE
echo "#PBS -o ./output.txt" >> $SUBMITFILE
echo  "" >> $SUBMITFILE

if [ $DECISION -eq 1 ]
then
    echo "#PBS -l nodes="$NUMBERNODES":ppn="$NUMBERCORES >> $SUBMITFILE
else
    echo "#PBS -l nodes=1:ppn=1" >> $SUBMITFILE
    echo "#PBS -t 0-"$ARRAYMAX >> $SUBMITFILE
    echo "ARRAYMAX="$ARRAYMAX >> $SUBMITFILE
    echo "export ARRAYMAX" >> $SUBMITFILE
fi

echo  "" >> $SUBMITFILE

if [ $DECISION -eq 1 ]
then
    echo -n "mpiexec --display-map -mca plm rsh -mca btl openib,self,sm -n " >> $SUBMITFILE
    echo  $(echo $NUMBERNODES " * " $NUMBERCORES | bc )" ./executable" >> $SUBMITFILE
    make MPI
else
    echo  " ./executable" >> $SUBMITFILE    
    make ARRAY
fi

echo  "" >> $SUBMITFILE




