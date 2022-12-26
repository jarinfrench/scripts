#! /bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=128
#SBATCH -t 24:00:00
#SBATCH -p normal_q
#SBATCH -A conc_alloys
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jarinf@vt.edu
#SBATCH -J UMoXe_gg_1

#SBATCH --array=516,533,577
##SBATCH --array=0-239%8
##SBATCH --array=240-479%8
##SBATCH --array=480-647%8

module reset
module load foss/2020a

LAMMPS_EXEC=/home/jarinf/LAMMPS/lammps/tinkercliffs_mpi_build_foss/lmp_29Oct2020
ffg=/home/jarinf/scripts/tinkercliffs_foss_2020a_executables/find_grains
cga=/home/jarinf/scripts/tinkercliffs_foss_2020a_executables/calculate_grain_area
potential_db=/home/jarinf/scripts/lattice_params.db
INFILE=/home/jarinf/projects/U/moly_xenon_effect/grain_growth/gamma/UMoXeVaried_MSD_gg.in
DATAFILES=(/home/jarinf/projects/U/{LAMMPS_U_100_75x75x5_pTernaryEAM_{20,30,45}.00degree_rcut0.70_l3.542_r99A_annealed.dat,LAMMPS_U_110_75x53x4_pTernaryEAM_{20,30,45}.00degree_l3.542_r99A_annealed.dat,LAMMPS_U_111_53x31x3_pTernaryEAM_{20,30,45}.00degree_rcut0.70_l3.542_r99A_annealed.dat})
Mo=(1 3 5)
Xe=(0.25 0.50 0.75)
T=({1050..1400..50})
timesteps=4000000 # I do not expect complete grain growth to occur e.g. in U5Mo0.75Xe, if any at all, at the lower temperatures. I am hoping to get some growth data at the higher temperatures, and expect it to take longer for these simulations, thus all will use the same timesteps value of 4000000.

T_idx=$((SLURM_ARRAY_TASK_ID%${#T[@]})) # iterate over 0-7 repeatedly
cXe_idx=$(($(echo "${SLURM_ARRAY_TASK_ID}/${#T[@]}" | bc) % ${#Xe[@]})) # iterate once every 8 values, looping about 0-2
cMo_idx=$(($(echo "${SLURM_ARRAY_TASK_ID}/(${#T[@]}*${#Xe[@]})" | bc) % ${#Mo[@]})) # iterate once every 3*8 values, looping about 0-2
df_idx=$(($(echo "${SLURM_ARRAY_TASK_ID}/(${#T[@]}*${#Xe[@]}*${#Mo[@]})" | bc) % ${#DATAFILES[@]})) # iterate once every 3*3*8 values, looping about 0-9

datafile=${DATAFILES[${df_idx}]}
axis=$(basename "${datafile}" | awk -F'_' '{print $3}')
cMo=${Mo[${cMo_idx}]}
cXe=${Xe[${cXe_idx}]}
mis_dir=$(basename "${datafile}" | awk -F'_' '{print $6}' | cut -c 3-5 --complement) # cut out the .00 from the file name, this gives us the misorientation directory (mis_dir)
Temp=${T[${T_idx}]}
dt=${timesteps}

datadir=/home/jarinf/projects/U/moly_xenon_effect/grain_growth/gamma/${axis}/ternary_eam/${cMo}at%Mo/${cXe}at%Xe/${mis_dir}/T${Temp}/large_r

for i in final_{1..3}; do
  (
  SEED=${RANDOM}
  newfile=$(basename "${datafile}" .dat | awk -v Mo="${cMo}" -v Xe="${cXe}" -F'U' '{print $1"U"Mo"Mo"Xe"Xe"$2".dat"}')
  vars_list="-var SEED ${SEED} -var file ${datafile} -var T ${Temp} -var timesteps ${dt} -var Mo ${cMo} -var Xe ${cXe} -var newfile ${newfile}"
  mkdir -p "${datadir}/dir_${i}"
  cd "${datadir}/dir_${i}" || exit
  mpirun ${LAMMPS_EXEC} ${vars_list} < ${INFILE} > "output_${Temp}.txt"
  #         num_files             misorientation                                num_types neigh_cutoff orientation_cutoff a0 structure
  echo "$(find "${datadir}/dir_${i}" -name  "*.dump" | wc -l) $(basename "${datafile}" | awk -F'_' '{print $6}' | cut -c 1-5) 3 1.207 0.0 3.542 bcc" > find_grains_input.txt
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
  ls -v ./*.dump >> find_grains_input.txt
  ${ffg} find_grains_input.txt -e 0
  echo "data.txt ${Temp} $(printf '%.3f\n' "$(echo "scale=3; ${cMo}/100" | bc)") $(awk 'NR==7 {print $2}' "${datafile}") 3.542 bcc" > grain_area_input.txt
  ${cga} grain_area_input.txt --potential-file ${potential_db} -p 3
  tar -cjvf snapshots.tar.bz2 ./*.dump && rm ./*.dump
  )
done

archive="$(dirname "$(dirname "${datadir}")")/T${Temp}"
echo "Archiving ${archive}"
tar -cjvf "/home/jarinf/archive_tmp/U${cMo}Mo${cXe}Xe_${axis}_${mis_dir}_T${Temp}_results.tar.bz2" "${archive}" && rm -rf "${archive}" && echo "Successfully archived U${cMo}Mo${cXe}Xe_${axis}_${mis_dir}_T${Temp}_results!"

