#! /bin/bash

# generates the pbs scripts with the accompanying LAMMPS input files

for i in {1..60}
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

mpirun -np $``PBS_NP /home/jarinf/LAMMPS/lammps-17Nov16/src/lmp_openmpigccfftw -var SEED \`bash -c 'echo $``RANDOM'\` < UO2_structure_minimize_${i}degree_rcut2.5.in > minimize_${i}degree_rcut2.5.txt

exit;" >> lmp_minimize_${i}degree_rcut2.5.pbs

  sed "24s[.*[read_data /home/jarinf/UO2_Circular_Grain/Atoms_removed/LAMMPS_UO2_SC_${i}degree_r100A_removed_rcut2.5.dat[" UO2_structure_minimization.in > UO2_structure_minimize_${i}degree_rcut2.5.in
  sed -i "136s[.*[dump atompos3 all custom 10000 dump3.pos.${i}degree_rcut2.5.*.dat id type q x y z c_pe_layer1[" UO2_structure_minimize_${i}degree_rcut2.5.in
done
