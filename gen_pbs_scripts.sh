#! /bin/bash

#read -p "What is the cutoff radius used? " rcut
read -p "What is the rotation axis? " axis
read -p "What is the grain radius? " radius

if [ $axis -eq 111 ]; then
  range=120
elif [ $axis -eq 110 ]; then
  range=180
elif [ $axis -eq 100 ]; then
  range=90
else
  range=360
fi

# NOTE: the directory location of the removed atoms files may change!

read -p "Generate scripts from a txt file? " usefile
case $usefile in
  y|Y)
    read -e -p "Please enter the filename containing the angles: " FN
    while read -r theta; do
      echo "#!/bin/bash

    # qsub script for Cascades

    # NOTE: You will need to edit the Walltime, Resource Request, Queue, and Module lines
    # to suit the requirements of your job. You will also, of course have to replace the example job

    #PBS -l nodes=1:ppn=32
    #PBS -l walltime=12:00:00
    #PBS -q normal_q
    #PBS -A FeCr_Bai
    #PBS -W group_list=cascades

    #PBS -M jarinf@vt.edu
    #PBS -m bea

    module purge
    module load gcc openmpi

    cd $``PBS_O_WORKDIR

    mpirun -np $``PBS_NP /home/jarinf/LAMMPS/lammps-17Nov16/src/lmp_openmpigccfftw -var SEED \`bash -c 'echo $``RANDOM'\` < UO2_structure_minimize_${axis}_${theta}degree.in > minimize_${axis}_${theta}degree.txt

    exit;" >> lmp_minimize_${axis}_${theta}degree.pbs

      sed "24s[.*[read_data /home/jarinf/UO2_Circular_Grain/Atoms_removed/${axis}Tilt/LAMMPS_UO2_SC_${axis}_${theta}degree_r${radius}A_removed.dat[" UO2_structure_minimization.in > UO2_structure_minimize_${axis}_${theta}degree.in
      sed -i "136s[.*[dump atompos3 all custom 10000 dump3.pos.${axis}.${theta}degree.*.dat id type q x y z c_pe_layer1[" UO2_structure_minimize_${axis}_${theta}degree.in
    done < "$FN"
    ;;
  n|N)
    echo "Generating the default files..."
    # generates the pbs scripts with the accompanying LAMMPS input files
    for i in $(seq 1 $range);
    do
      echo "#!/bin/bash

# qsub script for Cascades

# NOTE: You will need to edit the Walltime, Resource Request, Queue, and Module lines
# to suit the requirements of your job. You will also, of course have to replace the example job

#PBS -l nodes=1:ppn=32
#PBS -l walltime=8:00:00
#PBS -q normal_q
#PBS -A FeCr_Bai
#PBS -W group_list=cascades

#PBS -M jarinf@vt.edu
#PBS -m bea

module purge
module load gcc openmpi

cd $``PBS_O_WORKDIR

mpirun -np $``PBS_NP /home/jarinf/LAMMPS/lammps-17Nov16/src/lmp_openmpigccfftw -var SEED \`bash -c 'echo $``RANDOM'\` < UO2_structure_minimize_${axis}_${i}degree.in > minimize_${axis}_${i}degree.txt

exit;" >> lmp_minimize_${axis}_${i}degree.pbs

      sed "24s[.*[read_data /home/jarinf/UO2_Circular_Grain/Atoms_removed/${axis}Tilt/LAMMPS_UO2_SC_${axis}_${i}degree_r${radius}A_removed.dat[" UO2_structure_minimization.in > UO2_structure_minimize_${axis}_${i}degree.in
      sed -i "136s[.*[dump atompos3 all custom 10000 dump3.pos.${axis}.${i}degree.*.dat id type q x y z c_pe_layer1[" UO2_structure_minimize_${axis}_${i}degree.in
    done
    ;;
  *)
    echo "Unrecognized option.  Please enter y|Y or n|N."
esac
