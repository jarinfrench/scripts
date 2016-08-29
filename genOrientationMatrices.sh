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
  exit 1
fi

FN=$1 # This takes the first argument from the command line - this is assumed to be a filename of the format 100Tilt

echo "Determining the axis and type of misorientation..."
AXIS=`echo $FN | grep -o "1[01][01]"` # Pulls out the axis from the input file name
TYPE=`echo $FN | grep -o "T[a-z]\{3,4\}"` # Pulls out the misorientation type from the file name

# Because we want this automated, we are just going to assume a grain boundary normal
re='[Tt]wist'
if [[ $TYPE =~ $re ]]; then
  NORM=$AXIS
fi

re='[Tt]ilt'
if [[ $TYPE =~ $re ]]; then
  if [[ $AXIS =~ "100" ]]; then
    NORM="010"
  fi
  if [[ $AXIS =~ "110" ]]; then
    NORM="1-10"
  fi
  if [[ $AXIS =~ "111" ]]; then
    NORM="1-10"
  fi
fi

# Otherwise we can't run the program properly
if [ -z ${NORM+x} ]; then
  echo "Type not recognized"
  echo $NORM
  echo $TYPE
  exit 2
fi

echo "Reading the file..."
IFS=","
[ ! -f $FN ] && { echo "$FN file not found"; exit 99; }
echo "Running the command: ~/projects/scripts/orientation_matrix.py $AXIS <angle> $NORM -s -q"
while read -r angle en; do # read the file with comma separated values

  ~/projects/scripts/orientation_matrix.py $AXIS $angle $NORM -s -q
done < "$FN"

IFS=$OLDIFS
