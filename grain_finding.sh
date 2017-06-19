#! /bin/bash

# This script will parse the current directory and find all of the dump files.
# Note that this is looking specifically for files in the format *.dump!  This
# may change later.

targets=($(ls *.dump))
read -p "Please enter the misorientation angle for the dump files in this directory: " angle

for i in "${targets[@]}"
do
  find_grains $i $angle
done

# May change this to write the actual output to a file and parse it there.
