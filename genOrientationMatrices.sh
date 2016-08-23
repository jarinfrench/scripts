#! /bin/bash

# This script will generate the orientation matrices through python by looping
# through the CSV values given in the input files for the fitting function found
# in the fitting_5_12_16 folder.
# Argument(s):
#    $1       Should be a filename from the fitting_5_12_16 folder that specifies
#             the angles and relative energies for the 100, 110, and 111 symmetric
#             tilt and twist boundaries

if [ "$#" -ne 1 ]; then
  echo "Illegal number of parameters"
  exit
fi

FN=$1 # This takes the first argument from the command line - this is assumed to be a filename of the format 100Tilt

echo "Determining the axis and type of misorientation..."
AXIS=`echo $FN | grep -o "1[01][01]"` # Pulls out the axis from the input file name
TYPE=`echo $FN | grep -o "T[a-z]\{3,4\}"` # Pulls out the misorientation type from the file name

echo "Reading the file..."
IFS=","
[ ! -f $FN ] && { echo "$FN file not found"; exit 99; }
echo "Running the command: ~/projects/scripts/orientation_matrix.py $AXIS <angle> $TYPE -s -q"
while read -r angle en; do # read the file with comma separated values

  ~/projects/scripts/orientation_matrix.py $AXIS $angle $TYPE -s -q
done < "$FN"

IFS=$OLDIFS
