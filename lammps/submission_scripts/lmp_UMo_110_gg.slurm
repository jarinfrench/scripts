#! /bin/bash

#SBATCH --nodes=2
#SBATCH --ntasks=256
#SBATCH -t 120:00:00
#SBATCH -p normal_q
#SBATCH -A GB_Fuels
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jarinf@vt.edu
#SBATCH -J UMo_gg

module reset
module load foss/2020a

LAMMPS_EXEC=/home/jarinf/LAMMPS/lammps/tinkercliffs_mpi_build_foss/lmp
INFILE=/home/jarinf/projects/U/moly_effect/grain_growth/gamma/UMo_MSD_gg.in
DATAFILE=/home/jarinf/projects/U/moly_effect/grain_growth/gamma/110/ternary_eam/5at%/LAMMPS_U5Mo_110_pTernaryEAM_45degree_r100A_annealed.dat
T=<<T>>
timesteps=4000000

for i in 1; do
  (
  SEED=$(bash -c 'echo "${RANDOM}"')
  vars_list="-var SEED ${SEED} -var file ${DATAFILE} -var T ${T} -var timesteps ${timesteps}"
  mkdir dir_${i}
  cd dir_${i}
  mpirun ${LAMMPS_EXEC} ${vars_list} < ${INFILE} > output_${T}.txt
  )
done

