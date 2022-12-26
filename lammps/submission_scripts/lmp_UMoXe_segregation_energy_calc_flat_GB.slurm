#! /bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH -t 144:00:00

#SBATCH -p normal_q
#SBATCH -A GB_Fuels
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jarinf@vt.edu
#SBATCH -J UMo_segregation

module reset
module load foss/2020a

LAMMPS_EXEC=${HOME}/LAMMPS/lammps/tinkercliffs_mpi_build_foss/lmp_29Oct2020
INFILE=${HOME}/projects/U/segregation_energy_calc/flat/UMoXe_segregation_v3.in
track=${HOME}/scripts/tinkercliffs_foss_2020a_executables/track_atoms
d2d=${HOME}/scripts/tinkercliffs_foss_2020a_executables/dump2dat
DATAFILES=($(ls ${HOME}/projects/U/segregation_energy_calc/flat/U_*.lmp))
for file in "${DATAFILES[@]}"; do
  axis=$(basename ${file} | awk -F'_' '{print $2}')
  mis=$(basename ${file} | awk -F'_' '{print $3}' | cut -d . -f1)
  natoms=$(head -n3 ${file} | tail -n1 | awk '{print $1}')
  count=1
  for id in $(seq 1 ${natoms}); do
    imp=Mo
    vars="-var id ${id} -var datafile ${file} -var axis ${axis} -var mis ${mis} -var imp ${imp} -var run ${count}"
    settings="-in ${INFILE} -screen ${HOME}/projects/U/segregation_energy_calc/UMo_${axis}_${mis}_output.txt"
    mpirun ${LAMMPS_EXEC} ${vars} ${settings}
    imp=Xe
    vars="-var id ${id} -var datafile ${file} -var axis ${axis} -var mis ${mis} -var imp ${imp} -var run ${count}"
    settings="-in ${INFILE} -screen ${HOME}/projects/U/segregation_energy_calc/UXe_${axis}_${mis}_output.txt"
    mpirun ${LAMMPS_EXEC} ${vars} ${settings}
    ((++count))
  done
done
