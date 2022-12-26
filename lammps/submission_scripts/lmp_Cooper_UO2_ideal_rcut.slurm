#! /bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=128
#SBATCH -t 4:00:00
#SBATCH -p normal_q
#SBATCH -A GB_Fuels
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jarinf@vt.edu
#SBATCH -J uo2_eam_rcut
#SBATCH --array=0-5

module reset
module load foss/2020a

LAMMPS_EXEC=/home/jarinf/LAMMPS/lammps/tinkercliffs_mpi_build_foss/lmp_29Oct2020
INFILE=/home/jarinf/projects/uo2/ideal_rcut/Cooper/equilibrate.in
DATAFILE=(/home/jarinf/projects/uo2/ideal_rcut/Cooper/LAMMPS_UO2_{100_50x50x5,110_50x35x4,111_35x20x3}_pCooper.dat)
rotandrem=/home/jarinf/scripts/tinkercliffs_foss_2020a_executables/rotate_and_remove

if [[ "${SLURM_ARRAY_TASK_ID}" -lt 2 ]]; then
  f=${DATAFILE[0]}
elif [[ "${SLURM_ARRAY_TASK_ID}" -lt 4 ]]; then
  f=${DATAFILE[1]}
else
  f=${DATAFILE[2]}
fi

if [[ $((${SLURM_ARRAY_TASK_ID} % 2)) -eq 0 ]]; then
  angle=20.00
else
  angle=45.00
fi

for rcut in $(seq 0.05 0.05 0.95); do
  axis=$(echo ${f} | awk -F'_' '{print $4}')
  echo "${f} cooper_${axis}_rcut${rcut}_${angle}degree 0.555 ${angle} UO2 fluorite 5.454 ${rcut}" > rotate_and_remove_${SLURM_ARRAY_TASK_ID}.input
  ${rotandrem} rotate_and_remove_${SLURM_ARRAY_TASK_ID}.input --no-index --charge
  vars="-var data_file cooper_${axis}_rcut${rcut}_${angle}degree_removed.dat"
  mpirun ${LAMMPS_EXEC} ${vars} < ${INFILE} > output_${axis}_rcut${rcut}_${angle}degree_medium.txt
done

