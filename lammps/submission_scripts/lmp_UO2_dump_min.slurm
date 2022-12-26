#! /bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH -t 6:00:00
#SBATCH -p preemptable_q
#SBATCH -A personal
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jarinf@vt.edu
#SBATCH -J dump_min

module reset
module load foss/2020a

LAMMPS=/home/jarinf/LAMMPS/lammps/tinkercliffs_mpi_build_foss/lmp_29Oct2020
INFILE=/home/jarinf/projects/uo2/dump_min/dump_minimization.in
DATAFILE=/home/jarinf/projects/uo2/grain_growth/cylindrical/LAMMPS_UO2_111_35x20x3_pCooper_45.00degree_r100A_rcut0.45_annealed.dat

for i in LAMMPS_*; do
  vars="-var file ${i} -var SEED ${RANDOM} -in ${INFILE} -sc output.txt"
  mpirun ${LAMMPS} ${vars}
done
