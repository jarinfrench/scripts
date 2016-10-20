#! /bin/bash

# This script will generate the orientation matrices through python by looping
# through the CSV values given in the input files.
# Argument(s):
#    $1       Should be a filename that specifies the angles and relative
#             energies for the 100, 110, and 111 symmetric tilt and twist
#             boundaries

# Command-line argument counter that checks for the correct number of arguments.
# Does not check for correct syntax.
if [ "$#" -ne 1 ]; then
  echo "Illegal number of parameters"
  exit 1
fi

# This takes the first argument from the command line - this is assumed to be a
# filename of the format 100Tilt.
FN=$1

echo "Determining the axis..."
# Pulls out the axis from the input file name.  This uses regex syntax to find
# a series of numbers that match either 100, 110, or 111.  This also has an issue
# where it will find a match for 101, but as long as the files are named correctly
# it shouldn't be an issue.
AXIS=`echo $FN | grep -o "1[01][01]"`

echo "Reading the file..."
IFS="," # separation character is the comma
# Exit with error code 99 if unable to read the file
[ ! -f $FN ] && { echo "$FN file not found"; exit 99; }

# This makes the assumption that the file orientation_matrix.py has executable
# rights.
echo "Running the command: ~/projects/scripts/orientation_matrix.py $AXIS <angle> -s -q"
while read -r angle en; do # read the file with comma separated values

  ~/projects/scripts/orientation_matrix.py $AXIS $angle -s -q
done < "$FN" # the "$FN" is required if it's going to run properly!

IFS=$OLDIFS # go back to the old separation character based on the system value.
