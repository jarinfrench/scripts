#! /bin/bash

# This script parses all available lammps dump files in the current directory into
# a format that LAMMPS can read.  Dump files must be listed as such (i.e. dump.xxx)

targets=($(ls | grep ^dump*.*.dat))

read -p "Please enter the radius of the rotated grains: " radius

for i in "${targets[@]}"
do
  parse_lammps_dump $i $radius
done
