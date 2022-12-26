#! /bin/bash

#SBATCH --nodes=2
#SBATCH --ntasks=256
#SBATCH -t 3:00:00
#SBATCH -p normal_q
#SBATCH -A GB_Fuels
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jarinf@vt.edu
#SBATCH -J U_vac_gg

module reset
module load foss/2020a

LAMMPS_EXEC=/home/jarinf/LAMMPS/lammps/tinkercliffs_mpi_build_foss/lmp_29Oct2020
INFILE=/home/jarinf/projects/U/vacancy_effect/grain_growth/gamma/U_vac_MSD_gg.in
DATAFILE=/home/jarinf/projects/U/vacancy_effect/grain_growth/gamma/110/ternary_eam/0.25at%/45degree/LAMMPS_U0.25Xe_110_pTernaryEAM_45degree_r100A_annealed.dat
T=<>
timesteps=<>

for i in test_{1..3}; do
  (
  SEED=$(bash -c 'echo "${RANDOM}"')
  vars_list="-var SEED ${SEED} -var file ${DATAFILE} -var T ${T} -var timesteps ${timesteps}"
  mkdir dir_${i}
  cd dir_${i}
  mpirun ${LAMMPS_EXEC} ${vars_list} < ${INFILE} > output_${T}.txt
  )
done

