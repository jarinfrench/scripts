#! /bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=128
#SBATCH -t 144:00:00

#SBATCH -p normal_q
#SBATCH -A fecr_bai
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jarinf@vt.edu
#SBATCH -J UMo_segregation

module reset
module load foss/2020a

LAMMPS_EXEC=${HOME}/LAMMPS/lammps/tinkercliffs_mpi_build_foss/lmp_29Oct2020
INFILE=${HOME}/projects/U/segregation_energy_calc/UMoXe_segregation_v3.in
track=${HOME}/scripts/tinkercliffs_foss_2020a_executables/track_atoms
d2d=${HOME}/scripts/tinkercliffs_foss_2020a_executables/dump2dat
DATAFILES=($(ls ${HOME}/projects/U/LAMMPS_*_annealed.dat))
T=1100
timesteps=50000
for file in "${DATAFILES[@]}"; do
  axis=$(basename ${file} | awk -F'_' '{print $3}')
  mis=$(basename ${file} | awk -F'_' '{print $6}' | cut -c 3-5 --complement)
  min_file=$(basename ${file%.dat}_minimized.dat)
  ${track} ${file} --print-atom-ids > atom_ids.txt
  mpirun ${LAMMPS_EXEC} -var datafile ${file} -in minimize.in -sc ${HOME}/projects/U/segregation_energy_calc/output_$(basename ${file%.dat}_min.txt)
  ${d2d} ${file}_minimized.dump ${min_file} -e U
  count=1
  while read -u 9 id; do
    imp=Mo # 1 for Mo, 2 for Xe
    vars="-var id ${id} -var datafile ${min_file} -var SEED ${RANDOM} -var axis ${axis} -var mis ${mis} -var imp ${imp} -var run ${count} -var T ${T} -var timesteps ${timesteps}"
    settings="-in ${INFILE} -screen ${HOME}/projects/U/segregation_energy_calc/UMo_${axis}_${mis}_output.txt"
    mpirun ${LAMMPS_EXEC} ${vars} ${settings}
    imp=Xe
    vars="-var id ${id} -var datafile ${min_file} -var SEED ${RANDOM} -var axis ${axis} -var mis ${mis} -var imp ${imp} -var run ${count} -var T ${T} -var timesteps ${timesteps}"
    settings="-in ${INFILE} -screen ${HOME}/projects/U/segregation_energy_calc/UXe_${axis}_${mis}_output.txt"
    mpirun ${LAMMPS_EXEC} ${vars} ${settings}
    ((++count))
  done 9< atom_ids.txt
done
