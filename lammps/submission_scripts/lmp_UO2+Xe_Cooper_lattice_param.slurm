#! /bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH -t 6:00:00
#SBATCH -p normal_q
#SBATCH -A GB_Fuels
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jarinf@vt.edu
#SBATCH -J UXe_gg

module reset
module load foss/2020a

LAMMPS_EXEC=/home/jarinf/LAMMPS/lammps/tinkercliffs_mpi_build_foss/lmp_29Oct2020
INFILE=/home/jarinf/projects/uo2/potential_validation/Cooper/uo2_xe/lattice_param.in

for a in $(seq 5 0.001 6); do
  for c in $(seq 0.001 0.001 2.0); do
    SEED=${RANDOM}
    vars_list="-var a0 ${a} -var c ${c} -var SEED ${SEED} -log log_${a}A_$($c*100)Xe.lammps"
    mpirun ${LAMMPS_EXEC} ${vars_list} < ${INFILE} > output_${a}_${c}.txt
  done
done

