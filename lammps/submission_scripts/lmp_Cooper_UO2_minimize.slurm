#! /bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=128
#SBATCH -t 1:00:00
#SBATCH -p normal_q
#SBATCH -A GB_Fuels
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jarinf@vt.edu
#SBATCH -J structure_minimize 

#SBATCH --array=0-39

module reset
module load foss/2020a

LAMMPS_EXEC=/home/jarinf/LAMMPS/lammps/tinkercliffs_mpi_build_foss/lmp_29Oct2020
datafiles=($(ls ../LAMMPS_*_pCooper_*.dat))
datafile=${datafiles[$SLURM_ARRAY_TASK_ID]}
axis=$(basename "${datafile}" | awk -F'_' '{print $3}')
radius=$(basename "${datafile}" | awk -F'_' '{print $7}')
mis=$(basename "${datafile}" | awk -F'_' '{print $6}' | cut -c 3-5 --complement)
label="Cooper_UO2_${axis}_${mis}_${radius}"
vars="-var isCoreShell 0 -var str_label ${label} -var data_file ${datafile}"
settings="-log ${label}_log.lammps -in minimize.in -screen ${label}_output.txt -sf omp -pk omp 1"

mpirun ${LAMMPS_EXEC} ${vars} ${settings}
