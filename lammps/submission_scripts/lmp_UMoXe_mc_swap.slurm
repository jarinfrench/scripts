#! /bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH -t 24:00:00
#SBATCH -p normal_q
#SBATCH -A GB_Fuels
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jarinf@vt.edu
#SBATCH -J UMoXe_MC

module reset
module load foss/2020a

LAMMPS_EXEC=/home/jarinf/LAMMPS/lammps/tinkercliffs_mpi_build_foss/lmp_29Oct2020
infile=/home/jarinf/projects/U/monte_carlo/UMoXe_mc_swap.in
datafiles=($(ls /home/jarinf/projects/U/monte_carlo/*.dat))
Mo=(1 3 5 10 14 17 21)
Xe=(0.25 0.50 0.75 1.00)

for f in "${datafiles[@]}"; do
  axis=$(echo "$(basename ${f})" | awk -F'_' '{print $3}')
  mis=$(echo "$(basename ${f})" | awk -F'_' '{print $6}' | cut -c 3-5 --complement)
  for cMo in "${Mo[@]}"; do
    varlist="-var datafile ${f} -var c ${cMo} -var SEED ${RANDOM} -var imp 1 -var axis ${axis} -var mis ${mis}"
    mpirun ${LAMMPS_EXEC} ${varlist} < ${infile} > ${axis}_${mis}_Mo_output.txt
  done

  for cXe in "${Xe[@]}"; do
    varlist="-var datafile ${f} -var c ${cXe} -var SEED ${RANDOM} -var imp 2 -var axis ${axis} -var mis ${mis}"
    mpirun ${LAMMPS_EXEC} ${varlist} < ${infile} > ${axis}_${mis}_Xe_output.txt
  done

done
