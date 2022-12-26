#! /bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH -t 12:00:00
#SBATCH -p preemptable_q
#SBATCH -A personal
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jarinf@vt.edu
#SBATCH -J UO2Xe_eam_seg

#SBATCH --array=0-134

module reset
module load foss/2020a

LAMMPS_EXEC=/home/jarinf/LAMMPS/lammps/tinkercliffs_mpi_build_foss/lmp_29Oct2020
impgen=/home/jarinf/scripts/tinkercliffs_foss_2020a_executables/generate_impurities
DATAFILES=($(ls /home/jarinf/projects/uo2/impurity_effect/grain_growth/LAMMPS_*_minimized.dat))
INFILE=/home/jarinf/projects/uo2/impurity_effect/segregation_energy/segregation.in
rootdir=/home/jarinf/projects/uo2/impurity_effect/segregation_energy

fileno=$(($(echo "${SLURM_ARRAY_TASK_ID}/15" | bc)%9))
split_fileno=$(printf "%02d" $((SLURM_ARRAY_TASK_ID%15)))
potential="Cooper"
datafile=${DATAFILES[${fileno}]}
axis=$(basename ${datafile} | awk -F'_' '{print $3}')
mis=$(basename ${datafile} | awk -F'_' '{print $6}')
radius=$(basename ${datafile} | awk -F'_' '{print $7}')

while read -u 9 id; do
  echo "${datafile} ${rootdir}/tmp_${fileno}_${split_fileno}_${id}.dat UO2 ${id} 5.454" | tee ${rootdir}/tmp_${fileno}_${split_fileno}_impgen.txt >> ${rootdir}/impurity_generation_inputs.log
  ${impgen} ${rootdir}/tmp_${fileno}_${split_fileno}_impgen.txt -d s --no-molecule-removal
  vars="-var id ${id} -var fileno ${fileno} -var split_fileno ${split_fileno} -var potential ${potential} -var axis ${axis} -var mis ${mis} -var radius ${radius} -var rootdir ${rootdir} -var datafile ${rootdir}/tmp_${fileno}_${split_fileno}_${id}.dat -in ${INFILE} -sc ${rootdir}/output_${fileno}_${split_fileno}_${id}.txt"
  mpirun ${LAMMPS_EXEC} ${vars}
  rm ${rootdir}/tmp_${fileno}_${split_fileno}_${id}.dat
done 9< ${rootdir}/tmp_${fileno}_${split_fileno}

