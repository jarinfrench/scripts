#! /bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH -t 2:00:00
#SBATCH -p normal_q
#SBATCH -A GB_Fuels
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jarinf@vt.edu
#SBATCH -J gg_postprocess
#SBATCH --array=0-8

module reset
module load foss/2020a

ffg=/home/jarinf/scripts/tinkercliffs_foss_2020a_executables/find_grains
cga=/home/jarinf/scripts/tinkercliffs_foss_2020a_executables/calculate_grain_area
potential_db=/home/jarinf/scripts/lattice_params.db
ppdir=(/home/jarinf/projects/U/moly_effect/grain_growth/gamma/110/ternary_eam/3at%/45degree/T{1300,1350,1400}/large_r/dir_test_{4,5,6})

#if [[ ${SLURM_ARRAY_TASK_ID} -le 9 ]]; then # Mo
  p=3
  c=0.01
  DATAFILE=/home/jarinf/projects/U/moly_effect/grain_growth/gamma/110/ternary_eam/1at%/45degree/LAMMPS_U_110_pTernaryEAM_45degree_r100A_annealed.dat
#else # Xe
#  p=6
#  c=0.005
#  DATAFILE=/home/jarinf/projects/U/xenon_effect/grain_growth/gamma/110/ternary_eam/0.50at%/45degree/LAMMPS_U_110_pTernaryEAM_45degree_r100A_annealed.dat
#fi
T=$(echo ${ppdir[$SLURM_ARRAY_TASK_ID]} | awk -F'/' '{print $13}' | cut -c 2-)
cd ${ppdir[${SLURM_ARRAY_TASK_ID}]}
echo -e "$(ls *.dump | wc -l) 45.00 2 1.207 0.0 3.542 bcc\n0 0 1\n1 -1 0\n1 1 0" > find_grains_input.txt
ls -v *.dump >> find_grains_input.txt
${ffg} find_grains_input.txt -e 0
echo "data.txt ${T} ${c} $(awk 'NR==7 {print $2}' ${DATAFILE}) 3.542 bcc" > grain_area_input.txt
${cga} grain_area_input.txt --potential-file ${potential_db} -p ${p}
tar -cjvf snapshots.tar.bz2 *.dump

