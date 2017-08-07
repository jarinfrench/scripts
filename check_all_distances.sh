#! /bin/bash

for i in {1..60}
  do
    echo "Checking angle $i"
    ./check_distances Removed/LAMMPS_UO2_SC_${i}degree_r100A_removed_rcut2.dat 100
  done

for i in {1..60}
  do
    echo "Checking angle $i"
    ./check_distances Removed/LAMMPS_UO2_SC_${i}degree_r100A_removed_rcut2.5.dat 100
  done

for i in {1..60}
  do
    echo "Checking angle $i"
    ./check_distances Removed/LAMMPS_UO2_SC_${i}degree_r100A_removed_rcut3.dat 100
  done
