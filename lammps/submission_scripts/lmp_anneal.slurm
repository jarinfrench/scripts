#! /bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH -t 12:00:00
#SBATCH -p normal_q
#SBATCH -A GB_Fuels
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jarinf@vt.edu
#SBATCH -J anneal
#SBATCH --array=0-11

module reset
module load foss/2020a

LAMMPS_EXEC=/home/jarinf/LAMMPS/lammps/tinkercliffs_mpi_build_foss/lmp_29Oct2020
D2D=/home/jarinf/scripts/tinkercliffs_foss_2020a_executables/dump2dat
reg="_U(([0-9]*)Mo)?(O2)?(([0-9]*(\.[0-9]*)?)?Xe)?_"
anneals_list=($(ls /home/jarinf/projects/uo2/impurity_effect/grain_growth/*/LAMMPS_UO2_*_pBasak_*.dat) $(ls /home/jarinf/projects/uo2/impurity_effect/grain_growth/*/LAMMPS_UO2_*_pCooper_*.dat))
INPUT=/home/jarinf/projects/anneal/anneal.in
#idx=0

# generate the impurity systems
#for angle in {20,30,45}; do
#  for Mo in {1,3}; do
#    vars_list="-var Mo ${Mo} -var angle ${angle} -in subs.in"
#    mpirun ${LAMMPS_EXEC} ${vars_list} > output_${angle}degree_${Mo}Mo.txt
#  done
#done

file=${anneals_list[${SLURM_ARRAY_TASK_ID}]}
path=$(dirname ${file})
base=$(basename ${file} .dat)
if [[ ${SLURM_ARRAY_TASK_ID} -le 5 ]]; then
  p=1
else
  p=2
fi

vars_list="-var base ${base} -var p ${p} -var path ${path:-.} -var SEED ${RANDOM} -in ${INPUT} -l log_${base}.lammps"
mpirun ${LAMMPS_EXEC} ${vars_list} > output_${base}.txt
${D2D} ${base}_annealed.dump ${base}_annealed.dat -e $(echo ${base} | awk -F'_' '{print $2}')


# convert the dump files to LAMMPS input structures
#for i in *.dump; do
#  ${D2D} ${i} ${i/%dump/dat} -e $(echo ${i} | awk -F'_' '{print $2}')
#  anneals_list[${idx}]=${i/%dump/dat}
#  let idx++
#done

# anneal each converted file as above
#for i in ${anneals_list[@]}; do
#  SEED=$(bash -c 'echo "${RANDOM}"')
#  if [[ "${i}" =~ ${reg} ]]; then
#    Mo="${BASH_REMATCH[2]:-0}"
#    Xe="${BASH_REMATCH[4]:-0}"
#  fi
#  vars_list="-var SEED ${SEED} -var base ${i/%.dat/} -var Mo ${Mo} -var Xe ${Xe} -in anneal.in"
#  mpirun ${LAMMPS_EXEC} ${vars_list} > output_${i/%.dat/}.txt
#done

# convert the annealed dump files to LAMMPS input files
#for i in *_annealed.dump do
#  ${D2D} ${i} ${i/%dump/dat} -e $(echo ${i} | awk -F'_' '{print $2}')
#done

