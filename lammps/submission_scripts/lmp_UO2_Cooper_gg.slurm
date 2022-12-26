#! /bin/bash

#SBATCH --nodes=2
#SBATCH --ntasks=256
#SBATCH -t 48:00:00
#SBATCH -p normal_q
#SBATCH -A GB_Fuels
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jarinf@vt.edu
#SBATCH -J UO2Xe_eam_gg_1

#SBATCH --array=0-239%24
##SBATCH --array=240-479%24
##SBATCH --array=480-719%24
##SBATCH --array=720-839%24 # Basak only

module reset
module load foss/2020a

LAMMPS_EXEC=/home/jarinf/LAMMPS/lammps/tinkercliffs_mpi_build_foss/lmp_29Oct2020
ffg=/home/jarinf/scripts/tinkercliffs_foss_2020a_executables/find_grains
cga=/home/jarinf/scripts/tinkercliffs_foss_2020a_executables/calculate_grain_area
impgen=/home/jarinf/scripts/tinkercliffs_foss_2020a_executables/generate_impurities
potential_db=/home/jarinf/scripts/lattice_params.db
rootdir=/home/jarinf/projects/uo2/impurity_effect/grain_growth
INFILE=/home/jarinf/projects/uo2/impurity_effect/grain_growth/UO2XeVaried_MSD_Cooper_gg.in
potential="Cooper"
Xe=(0.25 0.50 0.75 1.00 1.25 1.50)

if [[ ${potential} == "Cooper" ]]; then
  #Cooper stuff
  DATAFILES=($(ls ${rootdir}/LAMMPS_*_pCooper_*annealed.dat | sort -k 1.87,1.88)) 
  p=18
  T=(3050 {3000..2000..100}) #Cooper
elif [[ "${potential}" == "Basak" ]]; then
  #Basak stuff
  DATAFILES=($(ls ${rootdir}/LAMMPS_*_pBasak_*annealed.dat | sort -k 1.86,1.87))
  p=19
  T=({3300..2000..100}) #Basak
else
  echo "Unrecognized potential: ${potential} (basename = ${datafile})"
  exit
fi

T_idx=$((SLURM_ARRAY_TASK_ID%${#T[@]}))
c_idx=$(($(echo "${SLURM_ARRAY_TASK_ID}/${#T[@]}" | bc) % ${#Xe[@]}))
df_idx=$(($(echo "${SLURM_ARRAY_TASK_ID}/(${#T[@]}*${#Xe[@]})" | bc) % ${#DATAFILES[@]}))
datafile=${DATAFILES[${df_idx}]}

timesteps=4000000 # Just a single timestep value, since we expect significant slow down with the presence of Xe

axis=$(basename "${datafile}" | awk -F'_' '{print $3}')
radius=$(basename "${datafile}" | awk -F'_' '{print $7}')
if [[ ${radius} == "r100A" ]]; then
  radius_dir="large_r"
elif [[ ${radius} == "r75A" ]]; then
  radius_dir="medium_r"
else
  echo "Unrecognized grain radius: ${radius} (basename = ${datafile})"
  exit
fi

  
c=${Xe[${c_idx}]}
mis_dir=$(basename "${datafile}" | awk -F'_' '{print $6}' | cut -c 3-5 --complement) # cut out the .00 from the file name, this gives us the misorientation directory (mis_dir)
Temp=${T[${T_idx}]}

datadir=${rootdir}/${axis}/${mis_dir}/${potential}/T${Temp}/${radius_dir}
echo "${datadir} being processed" >> /home/jarinf/projects/uo2/impurity_effect/grain_growth/slurm_log.txt 

#for i in final_{2..3}; do
for i in final_1; do
  (
  mkdir -p "${datadir}/dir_${i}"
  SEED=${RANDOM}
  impurityfile=$(basename "${datafile}" .dat | awk -v Xe="${c}" -F'UO2' '{print $1"UO2+"Xe"Xe"$2".dat"}')
  echo "${datafile} ${datadir}/${impurityfile} UO2 $(printf '%.4f\n' "$(echo "scale=4; ${c}/100" | bc)") 1 ${RANDOM} 5.454" | tee "input.tmp${SLURM_ARRAY_TASK_ID}" >> "${rootdir}/generate_impurities_input.log"
  ${impgen} "input.tmp${SLURM_ARRAY_TASK_ID}" -d s
  vars_list="-var SEED ${SEED} -var file ${datadir}/${impurityfile} -var T ${Temp} -var timesteps ${timesteps} -var cXe ${c} -var newfile $(echo "$(basename ${impurityfile} .dat)_annealed.dump")"
  cd "${datadir}/dir_${i}" || exit
  mpirun ${LAMMPS_EXEC} ${vars_list} < ${INFILE} > "output_${Temp}.txt"
  #         num_files             misorientation                                num_types neigh_cutoff orientation_cutoff a0 structure
  echo "$(find "${datadir}/dir_${i}" -name  "*.dump" | wc -l) $(basename "${datafile}" | awk -F'_' '{print $6}' | cut -c 1-5) 3 1.207 0.0 5.454 fcc" > find_grains_input.txt
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
  ${ffg} find_grains_input.txt -e 0 -i 2
  echo "data.txt ${Temp} $(printf '%.3f\n' "$(echo "scale=4; ${c}/100" | bc)") $(awk 'NR==7 {print $2}' "${datafile}") 5.454 fcc" > grain_area_input.txt
  ${cga} grain_area_input.txt --potential-file ${potential_db} -p ${p}
  tar -cjvf snapshots.tar.bz2 ./*.dump && rm ./*.dump
  )
done

archive="$(dirname "$(dirname "${datadir}")")/T${Temp}"
echo "Archiving ${archive}"
#tar -cjvf "/home/jarinf/archive_tmp/U${c}Xe_${axis}_${mis_dir}_T${Temp}_results.tar.bz2" "${archive}" && rm -rf "${archive}" && echo "Successfully archived U${c}Xe_${axis}_${mis_dir}_T${Temp}_results!"
tar -cjvf "/home/jarinf/archive_tmp/UO2+${c}Xe_${axis}_${mis_dir}_T${Temp}_${radius}_dir1_results.tar.bz2" "${archive}" && rm -rf "${archive}" && echo "Successfully archived UO2+${c}Xe_${axis}_${mis_dir}_T${Temp}_${radius}_dir1_results!"

