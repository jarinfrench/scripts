#! /bin/bash

# This script generates the pbs scripts and input files for minimizing UO2 structures.

#read -p "What is the cutoff radius used? " rcut
read -p "What is the rotation axis? " axis
read -p "What is the grain radius? " radius
read -p "What is the element? " element # get the desired element

if [[ "${element,,}" =~ ^(uo2|al|cu|au|ni)$ ]]; then
  element=${element,,}
  if [[ "$element" =~ ^uo2$ ]]; then
    element=${element^^}
  else
    element=${element^}
  fi
else
  echo "$element is not recognized.  Please enter either UO2, Cu, Al, Au, or Ni."
  exit 3
fi

read -p "Is this for INL or VT? " destination # get the destination

if [[ "${destination,,}" =~ ^(inl|vt)$ ]]; then
  destination=${destination,,}
else
  echo "${destination} not recognized.  Please enter either INL or VT."
  exit 4
fi

if [ $axis -eq 111 ]; then # three-fold symmetry with the 111 axis
  let range=120/5
elif [ $axis -eq 110 ]; then # two-fold symmetry with the 110 axis
  let range=180/5
elif [ $axis -eq 100 ]; then # four-fold symmetry with the 100 axis
  let range=90/5
else
  # This may be too high, but remains to be determined.  See the rotate_box.py
  # script in the /projects/school/Research/UO2 directory to check the
  # particular axis being rotated about.
  let range=360/5
fi

read -e -p "Please enter the unrotated LAMMPS file name: " oFN
# Generate the no_gb files
if [[ "$destination" =~ ^inl$ ]]; then # Using INL
  echo "#!/bin/bash

# qsub script for Falcon

#PBS -j oe
#PBS -l select=4:ncpus=36:mpiprocs=36
#PBS -l walltime=01:00:00
#PBS -P neams

#PBS -M jarinf@vt.edu
#PBS -m bea
#PBS -N ${element}_no_GB

source /etc/profile.d/modules.sh
module load OpenMPI

cd $``PBS_O_WORKDIR

mpirun /home/frenjari/projects/lammps/lammps-17Nov16/src/lmp_openmpigccfftw_falcon -var SEED \`bash -c 'echo $``RANDOM'\` < ${element}_structure_minimize_${axis}_no_GB.in > minimize_${axis}_no_GB.txt

exit;" > lmp_minimize_${axis}_no_GB.pbs

  sed "7s[.*[read_data ${oFN}[" ${element}_structure_minimization.in > ${element}_structure_minimize_${axis}_no_GB.in
  if [[ "${element}" =~ ^UO2$ ]]; then
    sed -i "78s[.*[dump atompos1 all custom 10000 dump1.pos.no_gb.dat id type x y z c_pe_layer1[" ${element}_structure_minimize_${axis}_no_GB.in
    sed -i "119s[.*[dump atompos3 all custom 10000 dump3.pos.${axis}.no_GB.*.dat id type x y z c_pe_layer1[" ${element}_structure_minimize_${axis}_no_GB.in
  else
    sed -i "50s[.*[dump atompos1 all custom 10000 dump1.pos.no_gb.dat id type x y z c_pe_layer1[" ${element}_structure_minimize_${axis}_no_GB.in
    sed -i "91s[.*[dump atompos3 all custom 10000 dump3.pos.${axis}.no_GB.*.dat id type x y z c_pe_layer1[" ${element}_structure_minimize_${axis}_no_GB.in
  fi

else # Using ARC
  echo "#!/bin/bash

# qsub script for Cascades

#PBS -l nodes=1:ppn=32
#PBS -l walltime=00:30:00
#PBS -q normal_q
#PBS -A FeCr_Bai
#PBS -W group_list=cascades

#PBS -M jarinf@vt.edu
#PBS -m bea
#PBS -N ${element}_no_GB

module purge
module load gcc openmpi

cd $``PBS_O_WORKDIR

mpirun -np $``PBS_NP /home/jarinf/LAMMPS/lammps-17Nov16/src/lmp_openmpigccfftw_cascades -var SEED \`bash -c 'echo $``RANDOM'\` < ${element}_structure_minimize_${axis}_no_GB.in > minimize_${axis}_no_GB.txt

exit;" >> lmp_minimize_${axis}_no_GB.pbs

  sed "7s[.*[read_data ${oFN}[" ${element}_structure_minimization.in > ${element}_structure_minimize_${axis}_no_GB.in
  if [[ "${element}" =~ ^UO2$ ]]; then
    sed -i "78s[.*[dump atompos1 all custom 10000 dump1.pos.no_gb.dat id type x y z c_pe_layer1[" ${element}_structure_minimize_${axis}_no_GB.in
    sed -i "119s[.*[dump atompos3 all custom 10000 dump3.pos.${axis}.no_GB.*.dat id type x y z c_pe_layer1[" ${element}_structure_minimize_${axis}_no_GB.in
  else
    sed -i "50s[.*[dump atompos1 all custom 10000 dump1.pos.no_gb.dat id type x y z c_pe_layer1[" ${element}_structure_minimize_${axis}_no_GB.in
    sed -i "91s[.*[dump atompos3 all custom 10000 dump3.pos.${axis}.no_GB.*.dat id type x y z c_pe_layer1[" ${element}_structure_minimize_${axis}_no_GB.in
  fi

fi
# Sometimes we have a file with specific angles we want to check.  This handles
# that case.
read -p "Generate scripts from a txt file? " usefile
case $usefile in
  y|Y)
    read -e -p "Please enter the filename containing the angles: " FN
    while read -r theta; do
      if [[ "$destination" =~ ^inl$ ]]; then
        echo "#!/bin/bash

# qsub script for Falcon

#PBS -j oe
#PBS -l select=4:ncpus=36:mpiprocs=36
#PBS -l walltime=01:00:00
#PBS -P neams

#PBS -M jarinf@vt.edu
#PBS -m bea
#PBS -N ${element}_${axis}_${theta}degree

source /etc/profile.d/modules.sh
module load OpenMPI

cd $``PBS_O_WORKDIR

mpirun /home/frenjari/projects/lammps/lammps-17Nov16/src/lmp_openmpigccfftw_falcon -var SEED \`bash -c 'echo $``RANDOM'\` < ${element}_structure_minimize_${axis}_${theta}degree.in > minimize_${axis}_${theta}degree.txt

exit;" >> lmp_minimize_${axis}_${theta}degree.pbs

        sed "7s[.*[read_data LAMMPS_${element^^}_SC_${axis}_${theta}degree_r${radius}A_removed.dat[" ${element}_structure_minimization.in > ${element}_structure_minimize_${axis}_${theta}degree.in
        if [[ "${element}" =~ ^UO2$ ]]; then
          sed -i "78s[.*[dump atompos1 all custom 10000 dump1.pos.${theta}degree.dat id type x y z c_pe_layer1[" ${element}_structure_minimize_${axis}_${theta}degree.in
          sed -i "119s[.*[dump atompos3 all custom 10000 dump3.pos.${axis}.${theta}degree.dat id type x y z c_pe_layer1[" ${element}_structure_minimize_${axis}_${theta}degree.in
        else
          sed -i "50s[.*[dump atompos1 all custom 10000 dump1.pos.${theta}degree.dat id type x y z c_pe_layer1[" ${element}_structure_minimize_${axis}_${theta}degree.in
          sed -i "91s[.*[dump atompos3 all custom 10000 dump3.pos.${axis}.${theta}degree.dat id type x y z c_pe_layer1[" ${element}_structure_minimize_${axis}_${theta}degree.in
        fi
      else
        echo "#!/bin/bash

# qsub script for Cascades

#PBS -l nodes=1:ppn=32
#PBS -l walltime=1:00:00
#PBS -q normal_q
#PBS -A FeCr_Bai
#PBS -W group_list=cascades

#PBS -M jarinf@vt.edu
#PBS -m bea
#PBS -N ${element}_${axis}_${theta}degree

module purge
module load gcc openmpi

cd $``PBS_O_WORKDIR

mpirun -np $``PBS_NP /home/jarinf/LAMMPS/lammps-17Nov16/src/lmp_openmpigccfftw_cascades -var SEED \`bash -c 'echo $``RANDOM'\` < ${element}_structure_minimize_${axis}_${theta}degree.in > minimize_${axis}_${theta}degree.txt

exit;" >> lmp_minimize_${axis}_${theta}degree.pbs

        sed "7s[.*[read_data LAMMPS_${element^^}_${axis}_${theta}degree_r${radius}A_removed.dat[" ${element}_structure_minimization.in > ${element}_structure_minimize_${axis}_${theta}degree.in
        if [[ "${element}" =~ ^UO2$ ]]; then
          sed -i "78s[.*[dump atompos1 all custom 10000 dump1.pos.${theta}degree.dat id type x y z c_pe_layer1[" ${element}_structure_minimize_${axis}_${theta}degree.in
          sed -i "119s[.*[dump atompos3 all custom 10000 dump3.pos.${axis}.${theta}degree.dat id type x y z c_pe_layer1[" ${element}_structure_minimize_${axis}_${theta}degree.in
        else
          sed -i "50s[.*[dump atompos1 all custom 10000 dump1.pos.${theta}degree.dat id type x y z c_pe_layer1[" ${element}_structure_minimize_${axis}_${theta}degree.in
          sed -i "91s[.*[dump atompos3 all custom 10000 dump3.pos.${axis}.${theta}degree.dat id type x y z c_pe_layer1[" ${element}_structure_minimize_${axis}_${theta}degree.in
        fi
      fi
      done < "$FN"
    ;;
  n|N)
    echo "Generating the default files..."
    # generates the pbs scripts with the accompanying LAMMPS input files
    for i in $(seq 1 $range);
    do
      if [[ "$destination" =~ ^inl$ ]]; then
        echo "#!/bin/bash

# qsub script for Falcon

#PBS -j oe
#PBS -l select=4:ncpus=36:mpiprocs=36
#PBS -l walltime=01:00:00
#PBS -P neams

#PBS -M jarinf@vt.edu
#PBS -m bea
#PBS -N ${element}_${axis}_$((${i}*5))degree

source /etc/profile.d/modules.sh
module load OpenMPI

cd $``PBS_O_WORKDIR

mpirun /home/frenjari/projects/lammps/lammps-17Nov16/src/lmp_openmpigccfftw_falcon -var SEED \`bash -c 'echo $``RANDOM'\` < ${element}_structure_minimize_${axis}_$((${i}*5))degree.in > minimize_${axis}_$((${i}*5))degree.txt

exit;" >> lmp_minimize_${axis}_$((${i}*5))degree.pbs

        sed "7s[.*[read_data LAMMPS_${element^^}_SC_${axis}_$((${i}*5))degree_r${radius}A_removed.dat[" ${element}_structure_minimization.in > ${element}_structure_minimize_${axis}_$((${i}*5))degree.in
        if [[ "${element}" =~ ^UO2$ ]]; then
          sed -i "78s[.*[dump atompos1 all custom 10000 dump1.pos.$((${i}*5))degree.dat id type x y z c_pe_layer1[" ${element}_structure_minimize_${axis}_$((${i}*5))degree.in
          sed -i "119s[.*[dump atompos3 all custom 10000 dump3.pos.${axis}.$((${i}*5))degree.*.dat id type x y z c_pe_layer1[" ${element}_structure_minimize_${axis}_$((${i}*5))degree.in
        else
          sed -i "50s[.*[dump atompos1 all custom 10000 dump1.pos.$((${i}*5))degree.dat id type x y z c_pe_layer1[" ${element}_structure_minimize_${axis}_$((${i}*5))degree.in
          sed -i "91s[.*[dump atompos3 all custom 10000 dump3.pos.${axis}.$((${i}*5))degree.*.dat id type x y z c_pe_layer1[" ${element}_structure_minimize_${axis}_$((${i}*5))degree.in
        fi
      else
        echo "#!/bin/bash

# qsub script for Cascades

#PBS -l nodes=1:ppn=32
#PBS -l walltime=8:00:00
#PBS -q normal_q
#PBS -A FeCr_Bai
#PBS -W group_list=cascades

#PBS -M jarinf@vt.edu
#PBS -m bea
#PBS -N ${element}_${axis}_$((${i}*5))degree

module purge
module load gcc openmpi

cd $``PBS_O_WORKDIR

mpirun -np $``PBS_NP /home/jarinf/LAMMPS/lammps-17Nov16/src/lmp_openmpigccfftw_cascades -var SEED \`bash -c 'echo $``RANDOM'\` < ${element}_structure_minimize_${axis}_$((${i}*5))degree.in > minimize_${axis}_$((${i}*5))degree.txt

exit;" >> lmp_minimize_${axis}_$((${i}*5))degree.pbs

        sed "7s[.*[read_data LAMMPS_${element^^}_${axis}_$((${i}*5))degree_r${radius}A_removed.dat[" ${element}_structure_minimization.in > ${element}_structure_minimize_${axis}_$((${i}*5))degree.in
        if [[ "${element}" =~ ^UO2$ ]]; then
          sed -i "78s[.*[dump atompos1 all custom 10000 dump1.pos.$((${i}*5))degree.dat id type x y z c_pe_layer1[" ${element}_structure_minimize_${axis}_$((${i}*5))degree.in
          sed -i "119s[.*[dump atompos3 all custom 10000 dump3.pos.${axis}.$((${i}*5))degree.*.dat id type x y z c_pe_layer1[" ${element}_structure_minimize_${axis}_$((${i}*5))degree.in
        else
          sed -i "50s[.*[dump atompos1 all custom 10000 dump1.pos.$((${i}*5))degree.dat id type x y z c_pe_layer1[" ${element}_structure_minimize_${axis}_$((${i}*5))degree.in
          sed -i "91s[.*[dump atompos3 all custom 10000 dump3.pos.${axis}.$((${i}*5))degree.*.dat id type x y z c_pe_layer1[" ${element}_structure_minimize_${axis}_$((${i}*5))degree.in
        fi
      fi
    done
    ;;
  *)
    echo "Unrecognized option.  Please enter y|Y or n|N."
esac

# move the files to the correct directory.  Assumes that you are in a directory
# that has as subdirectories the axis rotated about.
mv lmp_minimize_${axis}_* ${axis}/
mv ${element}_structure_minimize_${axis}* ${axis}/
