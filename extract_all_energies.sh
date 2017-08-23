#! /bin/bash

# This script utilizes the C++ script extract_energy.cpp and extracts the energy
# from all files with the format minimize_*.txt, and writes to the file
# specified by FN.  Note that this assumes that extract_energy is found in PATH!

# extract the .txt files
targets=($(ls -v | grep ^minimize_))
value=($(ls | grep -E "(^minimize[_0-3]*_no_GB)")) # This is the single grain value

read -p "Please enter the filename to be written to: " FN

# find the index in targets that has the value we want
for i in "${!targets[@]}"; do
  if [[ "${targets[$i]}" = "${value}" ]]; then
    j=$i;
  fi
done

# The energies will not be calculated correctly if the first value written to the
# file is not the single grain, so exit early.
# This line checks if the string is empty (null)
# See https://stackoverflow.com/questions/42111155/explanation-of-bash-if-command
if [ -z ${j+x} ]; then
  echo "Error finding initial energy configuration."
  exit 2
else
  extract_energy ${targets[$j]} $FN # gets the base value for a single grain
fi

# Extract the energy for each value.
for i in "${targets[@]}"
do
  if [[ $i = ${targets[$j]} ]]; then
    continue; # We don't want to double count the single grain energy.
  fi
  extract_energy $i $FN
done
