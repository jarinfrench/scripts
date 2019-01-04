#! /bin/bash
# This script will call on the script rotate_and_remove to generate UO2 structures.
# If the executable is not found, an error is thrown.  This script needs the
# original (unrotated) crystal, as well as the radius of the grain.

read -e -p "Please enter the filename of the original matrix: " FN

read -p "Please enter the radius of the rotated grain: " radius

read -p "Is this type cylinder (1) or sphere (2)? " boundary_type

FLAGS=''
if [[ "${boundary_type}" -eq 1 ]]; then
  b_type="cylinder"
elif [[ "${boundary_type}" -eq 2 ]]; then
  b_type="sphere"
  FLAGS='-s'
else
  echo "$boundary_type is not equal to 1 or 2.  Please enter 1 or 2 next time."
  exit 4
fi

read -p "Please enter the number of atom types: " ntypes

combinations=$(((($ntypes + 1) * $ntypes) / 2))

read -p "Please enter the cutoff radii separated by a space: " cutoff

countElems() { echo $#;}
numElems=`countElems $cutoff`
if [[ $numElems -ne $combinations ]]; then
  echo "Please enter $combinations radii for this system."
  exit 5
fi

echo $FN $radius $ntypes $cutoff > rotate_input.txt


# based on the filename, the axis is determined.  Can only handle 100 to 135 at
# this point.  Further modifications may be necessary to handle larger axes.
# Currently assumes (as of 1 March 2018) that the axis is immediately before a
# period, and that no other set of three numbers immediately precedes a period.
axis=`echo $FN | grep -o "[0-9][0-9][0-9]\." | cut -d. -f1`

# Utilize the symmetry of the axis.
if [ ${axis} -eq 111 ]; then
  let range=120/5
elif [ ${axis} -eq 110 ]; then
  let range=180/5
elif [ ${axis} -eq 100 ]; then
  let range=360/5
else
  echo "Not using a high-symmetry axis: setting range = 360 degrees"
  let range=360/5 # This may be too high, but generally not.
fi

# If we have a text file with the angles we are interested in checking, specify
# that here.
read -p "Are you reading angles from a txt file? " angles

case ${angles} in
  y|Y)
    read -e -p "Please enter the filename of the angles txt file: " FN2
    while read -r theta; do
      rotate_and_remove rotate_input.txt ${theta} ${FLAGS}
    done < "${FN2}"
    for i in $(seq 1 ${range});
    do
      rotate_and_remove rotate_input.txt $(($i*5)) ${FLAGS}
    done
    ;;
  n|N)
    echo "Generating default layouts..."
    for i in $(seq 1 ${range});
    do
      rotate_and_remove rotate_input.txt $(($i*5)) ${FLAGS}
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

# Make the directories
mkdir -p Removed/${axis} Marked/${axis} Rotated/${axis}

# Assumes that these directories already exist.
mv *_removed*.dat Removed/${axis}
mv *_marked*.dat Marked/${axis}
mv *_rotated*.dat Rotated/${axis}
