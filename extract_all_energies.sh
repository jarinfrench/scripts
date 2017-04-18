#! /bin/bash

# This script utilizes the C++ script extract_energy.cpp and extracts the energy
# from all files with the format minimize_*.txt, and writes to the file
# specified by FN.

# extract the .txt files
targets=($(ls | grep ^minimize_))
value=($(ls | grep -E "(^minimize_[0-3]*_no_GB)")) # This is the single grain value

read -p "Please enter the filename to be written to: " FN

# find the index in targets that has the value we want
for i in "${!targets[@]}"; do
  if [[ "${targets[$i]}" = "${value}" ]]; then
    j=$i;
  fi
done

# The energies will not be calculated correctly if the first value written to the
# file is not the single grain, so exit early.
if [ -z ${j+x} ]; then
  echo "Error finding initial energy configuration."
  exit 2
else
  #echo "Running command: ./extract_energy ${targets[$j]} $FN"
  ./extract_energy ${targets[$j]} $FN # gets the base value for a single grain
fi

# Extract the energy for each value.
for i in "${targets[@]}"
do
  if [[ $i = ${targets[$j]} ]]; then
    continue; # We don't want to double count the single grain energy.
  fi
  #echo "Running command: ./extract_energy $i $FN"
  ./extract_energy $i $FN
done
