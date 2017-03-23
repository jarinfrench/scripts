#! /bin/bash
# This script will first generate the rotated grains in increments of 1 degree
# based on the input file given.  If no input file is given, a prompt is provided.
# After rotating the grains, this script will call the remove_extras script.
# Errors will be thrown if the executables rotate_grain and remove_extras are not
# found.  The resulting files will be saved in the Removed and Rotated directories

read -e -p "Please enter the filename of the original UO2 matrix: " FN

read -p "Please enter the radius of the rotated grain: " radius

# Make sure we're in the correct directory
cd ~/projects/school/Research/UO2

for i in {1..60};
do
  if [ $i -eq 1 ]; then
    echo "Rotating $i degree"
  else
    echo "Rotating $i degrees"
  fi
  ./rotate_and_remove $FN $radius $i
done

# Now move all of the files with the removed atoms to the Removed directory,
# and the atoms with just the rotated atoms to the rotated directory
# Also move the file with marked atoms to the Marked directory (useful for
# visualization in TecPlot)
mv *_marked*.dat Marked/
mv *_removed*.dat Removed/
mv *_rotated*.dat Rotated/
