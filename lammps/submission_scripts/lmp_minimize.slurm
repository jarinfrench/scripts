#! /bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=128
#SBATCH -t 1:00:00
#SBATCH -p normal_q
#SBATCH -A GB_Fuels
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jarinf@vt.edu
#SBATCH -J structure_minimize 

module reset
module load foss/2020a

LAMMPS_EXEC=/home/jarinf/LAMMPS/lammps/tinkercliffs_mpi_build_foss/lmp
label="Cooper_111_sigma7"
vars="-var isCoreShell 0 -var str_label ${label} -var data_file LAMMPS_CeO2_111_35x20x3_pCooper_38.20degree_r97A_l5.3972_core.lmp"
settings="-log ${label}_log.lammps -in minimize.in -screen ${label}_output.txt -sf omp -pk omp 1"

mpirun ${LAMMPS_EXEC} ${vars} ${settings}
