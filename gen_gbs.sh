#! /bin/bash
# This script will call on the script rotate_and_remove to generate UO2 structures.
# If the executable is not found, an error is thrown.  This script needs the
# original (unrotated) crystal, as well as the radius of the grain.

read -e -p "Please enter the filename of the original matrix: " FN

read -p "Please enter the radius of the rotated grain: " radius

read -p "Is this type cylinder (1) or sphere (2)? " boundary_type

if [[ "${boundary_type}" -eq 1 ]]; then
  b_type="cylinder"
elif [[ "${boundary_type}" -eq 2 ]]; then
  b_type="sphere"
else
  echo "$boundary_type is not equal to 1 or 2.  Please enter 1 or 2 next time."
  exit 4
fi

# based on the filename, the axis is determined.  Can only handle 100 to 133 at
# this point.  Further modifications may be necessary to handle larger axes.
axis=`echo $FN | grep -o "1[0-3][0-5]"`

# Utilize the symmetry of the axis.
if [ $axis -eq 111 ]; then
  let range=120/5
elif [ $axis -eq 110 ]; then
  let range=180/5
elif [ $axis -eq 100 ]; then
  let range=90/5
else
  let range=360/5 # This may be too high, but generally not.
fi

# If we have a text file with the angles we are interested in checking, specify
# that here.
read -p "Are you reading angles from a txt file? " angles

case $angles in
  y|Y)
    read -e -p "Please enter the filename of the angles txt file: " FN2
    while read -r theta; do
      echo "Rotating $theta degrees"
      rotate_and_remove $FN $radius $theta $b_type
    done < "$FN2"
    ;;
  n|N)
    echo "Generating default layouts..."
    for i in $(seq 1 $range);
    do
      echo "Rotating $(($i*5)) degrees"
      rotate_and_remove $FN $radius $(($i*5)) $b_type
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

# Assumes that these directories already exist.
mv *_marked*.dat Marked/
mv *_removed*.dat Removed/
mv *_rotated*.dat Rotated/
