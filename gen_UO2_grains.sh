#! /bin/bash
# This script will first generate the rotated grains in increments of 1 degree
# based on the input file given.  If no input file is given, a prompt is provided.
# After rotating the grains, this script will call the remove_extras script.
# Errors will be thrown if the executables rotate_grain and remove_extras are not
# found.  The resulting files will be saved in the Removed and Rotated directories

read -e -p "Please enter the filename of the original UO2 matrix: " FN

read -p "Please enter the radius of the rotated grain: " radius

read -p "Are you reading angles from a txt file? " angles

case $angles in
  y|Y)
    read -e -p "Please enter the filename of the angles txt file: " FN2
    while read -r theta; do
      echo "Rotating $theta degrees"
      echo "  r_cut = 2.0"
      ./rotate_and_remove_rcut2 $FN $radius $theta
      echo "  r_cut = 2.5"
      ./rotate_and_remove_rcut2.5 $FN $radius $theta
      echo "  r_cut = 3.0"
      ./rotate_and_remove_rcut3 $FN $radius $theta
    done < "$FN2"
    ;;
  n|N)
    echo "Generating default layouts..."
    for i in {1..60};
    do
      if [ $i -eq 1 ]; then
        echo "Rotating $i degree"
      else
        echo "Rotating $i degrees"
      fi
      echo "  r_cut = 2.0"
      ./rotate_and_remove_rcut2 $FN $radius $i
      echo "  r_cut = 2.5"
      ./rotate_and_remove_rcut2.5 $FN $radius $i
      echo "  r_cut = 3.0"
      ./rotate_and_remove_rcut3 $FN $radius $i
    done
    ;;
  *)
    echo "Unrecognized option.  Please enter y|Y or n|N. Generating default layouts."
    ;;
esac

# Now move all of the files with the removed atoms to the Removed directory,
# and the atoms with just the rotated atoms to the rotated directory
# Also move the file with marked atoms to the Marked directory (useful for
# visualization in TecPlot)
mv *_marked*.dat Marked/
mv *_removed*.dat Removed/
mv *_rotated*.dat Rotated/
