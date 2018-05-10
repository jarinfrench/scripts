#! /bin/bash

cg1_targets=`echo $(ls -v CG1*.0.dump)`
cg2_targets=`echo $(ls -v CG2*)`

a=0
for i in $cg1_targets; do
  new_file=`echo $i | sed 's|CG1.\([0-9][0-9]*\(\.[0-9][0-9]\)\?degree\).0.dump|\1|'`"_tracked.dat"
  set -- $cg2_targets; shift $a; j=$1
  # now perform the calculation
  echo "Running command find_new_positions $i $j $new_file"
  find_new_positions $i $j $new_file
  a=`expr $a + 1`
done
