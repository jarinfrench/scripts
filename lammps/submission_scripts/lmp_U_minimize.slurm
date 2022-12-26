#! /bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=128
#SBATCH -t 4:00:00
#SBATCH -p normal_q
#SBATCH -A GB_Fuels
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jarinf@vt.edu
#SBATCH -J UMo_structure_minimize 

module reset
module load foss/2020a

LAMMPS_EXEC=/home/jarinf/LAMMPS/lammps/tinkercliffs_mpi_build_foss/lmp_29Oct2020

for i in 1 3 5 10 14 17 21; do
  for j in {20,30,45}degree; do
    label="U${i}Mo_${j}"
    vars="-var isCoreShell 0 -var str_label ${label} -var data_file LAMMPS_U${i}Mo_110_pTernaryEAM_${j}_r100A_annealed.dat"
    settings="-log ${label}_log.lammps -in minimize.in -screen ${label}_output.txt -sf omp -pk omp 1"
    mpirun ${LAMMPS_EXEC} ${vars} ${settings}
  done
done
