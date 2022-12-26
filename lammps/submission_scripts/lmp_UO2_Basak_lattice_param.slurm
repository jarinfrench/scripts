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
INFILE=/home/jarinf/projects/uo2/potential_validation/Basak/lattice_param/lattice_param.in

for a in $(seq 5 0.001 6); do
  vars_list="-var a0 ${a}"
  mpirun ${LAMMPS_EXEC} ${vars_list} < ${INFILE} > output_${a}.txt
done

