#! /bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=128
#SBATCH -t 1:30:00
#SBATCH -p normal_q
#SBATCH -A GB_Fuels
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jarinf@vt.edu
#SBATCH -J structure_minimize 

#SBATCH --array=0-1

module reset
module load foss/2020a

LAMMPS_EXEC=/home/jarinf/LAMMPS/lammps/tinkercliffs_mpi_build_foss/lmp_29Oct2020
d2d=/home/jarinf/scripts/tinkercliffs_foss_2020a_executables/dump2dat
datafiles=($(ls LAMMPS_UO2_110_*_pBasak_*r75A_rcut0.[0-9][0-9].dat))
datafile=${datafiles[$SLURM_ARRAY_TASK_ID]}
axis=$(basename "${datafile}" | awk -F'_' '{print $3}')
radius=$(basename "${datafile}" | awk -F'_' '{print $7}')
mis=$(basename "${datafile}" | awk -F'_' '{print $6}' | cut -c 3-5 --complement)
label="Basak_UO2_${axis}_${mis}_${radius}"
vars="-var isCoreShell 0 -var str_label ${label} -var data_file ${datafile}"
settings="-log ${label}_log.lammps -in minimize.in -screen ${label}_output.txt -sf omp -pk omp 1"

mpirun ${LAMMPS_EXEC} ${vars} ${settings}
tar -cjvf ${label}_minimization_results.tar.bz2 ${label}_*
${d2d} ${label}_final_min.dump ${datafile//.dat}_annealed.dat -e UO2
rm ${label}_*.dump
