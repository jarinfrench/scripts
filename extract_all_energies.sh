#! /bin/bash

# extract the .txt files
targets=($(ls | grep ^minimize_))
value='minimize_no_GB.txt' # holds the default values

# find the index in targets that has the value we want
for i in "${!targets[@]}"; do
  if [[ "${targets[$i]}" = "${value}" ]]; then
    j=$i;
  fi
done
#echo "Running command: ./extract_energy ${targets[$j]}"
./extract_energy ${targets[$j]} # gets the base value for a single grain

# Extract the energy for each value.
for i in "${targets[@]}"
do
  if [[ $i = ${targets[$j]} ]]; then
    continue;
  fi
  #echo "Running command: ./extract_energy $i"
  ./extract_energy $i
done
