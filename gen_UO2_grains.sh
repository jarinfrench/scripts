#! /bin/bash
# This script will first generate the rotated grains in increments of 1 degree
# based on the input file given.  If no input file is given, a prompt is provided.
# After rotating the grains, this script will call the remove_extras script.
# Errors will be thrown if the executables rotate_grain and remove_extras are not
# found.  The resulting files will be saved in the Removed and Rotated directories

read -e -p "Please enter the filename of the original UO2 matrix: " FN

read -p "Please enter the radius of the rotated grain: " radius

axis=`echo $FN | grep -o "1[0-3][0-3]"`

if [ $axis -eq 111 ]; then
  range=120
elif [ $axis -eq 110 ]; then
  range=180
elif [ $axis -eq 100 ]; then
  range=90
else
  #echo "Invalid rotation axis."
  #exit -1
  range=360 # This may not be how high we need to go - TODO: determine this!
fi

read -p "Are you reading angles from a txt file? " angles

case $angles in
  y|Y)
    read -e -p "Please enter the filename of the angles txt file: " FN2
    while read -r theta; do
      echo "Rotating $theta degrees"
      # note that this version of the executable creates a filename with 2
      # decimal points of precision for theta. (i.e. 90.00degree instead of 90degree)
      ./rotate_and_remove_v2 $FN $radius $theta
    done < "$FN2"
    ;;
  n|N)
    echo "Generating default layouts..."
    for i in $(seq 1 $range);
    do
      if [ $i -eq 1 ]; then
        echo "Rotating $i degree"
      else
        echo "Rotating $i degrees"
      fi
      ./rotate_and_remove $FN $radius $i
    done
    ;;
  *)
    echo "Unrecognized option.  Please enter y|Y or n|N."
    ;;
esac

# Now move all of the files with the removed atoms to the Removed directory,
# and the atoms with just the rotated atoms to the rotated directory
# Also move the file with marked atoms to the Marked directory (useful for
# visualization in TecPlot)
mv *_marked*.dat Marked/
mv *_removed*.dat Removed/
mv *_rotated*.dat Rotated/
