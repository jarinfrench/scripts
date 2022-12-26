#! /bin/bash

#SBATCH --nodes=2
#SBATCH --ntasks=256
#SBATCH -t 24:00:00
#SBATCH -p normal_q
#SBATCH -A scc_alloys
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jarinf@vt.edu
#SBATCH -J UMo_gg

#SBATCH --array=0-239%8
##SBATCH --array=240-479%8
##SBATCH --array=480-503%8

module reset
module load foss/2020a

LAMMPS_EXEC=/home/jarinf/LAMMPS/lammps/tinkercliffs_mpi_build_foss/lmp_29Oct2020
ffg=/home/jarinf/scripts/tinkercliffs_foss_2020a_executables/find_grains
cga=/home/jarinf/scripts/tinkercliffs_foss_2020a_executables/calculate_grain_area
potential_db=/home/jarinf/scripts/lattice_params.db
INFILE=/home/jarinf/projects/U/moly_effect/grain_growth/gamma/UMoVaried_MSD_gg.in
DATAFILES=(/home/jarinf/projects/U/{LAMMPS_U_100_75x75x5_pTernaryEAM_{20,30,45}.00degree_rcut0.70_l3.542_r99A_annealed.dat,LAMMPS_U_110_75x53x4_pTernaryEAM_{20,30,45}.00degree_l3.542_r99A_annealed.dat,LAMMPS_U_111_53x31x3_pTernaryEAM_{20,30,45}.00degree_rcut0.70_l3.542_r99A_annealed.dat})
Mo=(1 3 5 10 14 17 21)
T=({1050..1400..50})
timesteps=({4000000..1200000..400000}) # works for 21at%, so this will be fine for all others, at the expense of wasted simulation time

T_idx=$((SLURM_ARRAY_TASK_ID%${#T[@]}))
c_idx=$(($(echo "${SLURM_ARRAY_TASK_ID}/${#T[@]}" | bc) % ${#Mo[@]}))
df_idx=$(($(echo "${SLURM_ARRAY_TASK_ID}/(${#T[@]}*${#Mo[@]})" | bc) % ${#DATAFILES[@]}))
dt_idx=$((SLURM_ARRAY_TASK_ID%${#timesteps[@]}))

datafile=${DATAFILES[${df_idx}]}
axis=$(basename "${datafile}" | awk -F'_' '{print $3}')
c=${Mo[${c_idx}]}
mis_dir=$(basename "${datafile}" | awk -F'_' '{print $6}' | cut -c 3-5 --complement) # cut out the .00 from the file name, this gives us the misorientation directory (mis_dir)
Temp=${T[${T_idx}]}
dt=${timesteps[${dt_idx}]}

datadir=/home/jarinf/projects/U/moly_effect/grain_growth/gamma/${axis}/ternary_eam/${c}at%/${mis_dir}/T${Temp}/large_r
echo "${datadir} being processed" >> /home/jarinf/projects/U/slurm_log.txt

for i in final_{1..3}; do
  (
  SEED=${RANDOM}
  newfile=$(basename "${datafile}" .dat | awk -v Mo="${c}" -F'U' '{print $1"U"Mo"Mo"$2".dat"}')
  vars_list="-var SEED ${SEED} -var file ${datafile} -var T ${Temp} -var timesteps ${dt} -var Mo ${c} -var newfile ${newfile}"
  mkdir -p "${datadir}/dir_${i}"
  cd "${datadir}/dir_${i}" || exit
  mpirun ${LAMMPS_EXEC} ${vars_list} < ${INFILE} > "output_${Temp}.txt"
  #         num_files             misorientation                                num_types neigh_cutoff orientation_cutoff a0 structure
  echo "$(find "${datadir}/dir_${i}" -name  "*.dump" | wc -l) $(basename "${datafile}" | awk -F'_' '{print $6}' | cut -c 1-5) 2 1.207 0.0 3.542 bcc" > find_grains_input.txt
  case ${axis} in
    100)
      echo -e "1 0 0\n0 1 0\n0 0 1" >> find_grains_input.txt
      ;;
    110)
      echo -e "0 0 1\n1 -1 0\n1 1 0" >> find_grains_input.txt
      ;;
    111)
      echo -e "1 -1 0\n1 1 -2\n1 1 1" >> find_grains_input.txt
      ;;
    *)
      echo "ERROR: Unrecognized axis ${axis}"
      exit
      ;;
  esac
  ls -v ${datadir}/dir_${i}/*.dump >> find_grains_input.txt
  ${ffg} find_grains_input.txt -e 0
  echo "data.txt ${Temp} $(printf '%.3f\n' "$(echo "scale=3; ${c}/100" | bc)") $(awk 'NR==7 {print $2}' "${datafile}") 3.542 bcc" > grain_area_input.txt
  ${cga} grain_area_input.txt --potential-file ${potential_db} -p 3
  tar -cjvf snapshots.tar.bz2 ./*.dump && rm ./*.dump
  )
done

archive="$(dirname "$(dirname "${datadir}")")/T${Temp}"
echo "Archiving ${archive}"
tar -cjvf "/home/jarinf/archive_tmp/U${c}Mo_${axis}_${mis_dir}_T${Temp}_results.tar.bz2" "${archive}" && rm -rf "${archive}" && echo "Successfully archived U${c}Mo_${axis}_${mis_dir}_T${Temp}_results!"

