#! /bin/bash

#SBATCH --nodes=2
#SBATCH --ntasks=256
#SBATCH -t 36:00:00
#SBATCH -p preemptable_q
#SBATCH -A personal
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jarinf@vt.edu
#SBATCH -J UO2_Basak_gg

#SBATCH --array=0-104

module reset
module load foss/2020a

LAMMPS_EXEC=/home/jarinf/LAMMPS/lammps/tinkercliffs_mpi_build_foss/lmp_29Oct2020
ffg=/home/jarinf/scripts/tinkercliffs_foss_2020a_executables/find_grains
cga=/home/jarinf/scripts/tinkercliffs_foss_2020a_executables/calculate_grain_area
potential_db=/home/jarinf/scripts/lattice_params.db
INFILE=/home/jarinf/projects/uo2/grain_growth/cylindrical/UO2_MSD_Basak_gg.in
datafile_dir=/home/jarinf/projects/uo2/grain_growth/cylindrical
DATAFILES=(${datafile_dir}/LAMMPS_UO2_100_50x50x5_pBasak_20.00degree_r75A_rcut0.35_annealed.dat
${datafile_dir}/LAMMPS_UO2_110_50x35x4_pBasak_20.00degree_r75A_rcut0.50.dat
${datafile_dir}/LAMMPS_UO2_110_50x35x4_pBasak_45.00degree_r75A_rcut0.30.dat)

dirr=(1 2 3)

# note that T_idx is the index for both the temperature and the timesteps
if [[ ${SLURM_ARRAY_TASK_ID} -lt 42 ]]; then
  T=({2000..3300..100})
  timesteps=(4000000 4000000 4000000 3500000 3000000 2000000 1500000 1500000 1500000 1000000 1000000 800000 800000 800000)
  datafile=${DATAFILES[0]}
  T_idx=$(($(echo "${SLURM_ARRAY_TASK_ID}/(${#dirr[@]})" | bc) % ${#T[@]}))
elif [[ ${SLURM_ARRAY_TASK_ID} -lt 78 ]]; then
  T=({2200..3300..100})
  timesteps=(4000000 4000000 4000000 4000000 4000000 3500000 3500000 3500000 3500000 3500000 3500000 3500000)
  datafile=${DATAFILES[1]}
  T_idx=$(($(echo "(${SLURM_ARRAY_TASK_ID}-42)/(${#dirr[@]})" | bc) % ${#T[@]}))
else
  T=({2500..3300..100})
  timesteps=(4000000 4000000 4000000 4000000 4000000 3000000 2000000 1500000 1000000)
  datafile=${DATAFILES[2]}
  T_idx=$(($(echo "(${SLURM_ARRAY_TASK_ID}-78)/(${#dirr[@]})" | bc) % ${#T[@]}))
fi

axis=$(basename "${datafile}" | awk -F'_' '{print $3}')

radius=$(basename "${datafile}" | awk -F'_' '{print $7}')
radius_dir="medium_r"

mis_dir=$(basename "${datafile}" | awk -F'_' '{print $6}' | cut -c 3-5 --complement) # cut out the .00 from the file name, this gives us the misorientation directory (mis_dir)
Temp=${T[${T_idx}]}
dt=${timesteps[${T_idx}]}

datadir=$(dirname ${datafile})/${axis}/${mis_dir}/Basak/T${Temp}/${radius_dir}
echo "${datadir} being processed" >> /home/jarinf/projects/uo2/grain_growth/cylindrical/slurm_log.txt 
ii=$((SLURM_ARRAY_TASK_ID % ${#dirr[@]}))
i=${dirr[${ii}]}

#for i in final_{1..3}; do
  (
  SEED=${RANDOM}
  vars_list="-var SEED ${SEED} -var file ${datafile} -var T ${Temp} -var timesteps ${dt} -in ${INFILE} -sc ${datadir}/dir_final_${i}/output_${Temp}.txt"
  mkdir -p "${datadir}/dir_final_${i}"
  cd "${datadir}/dir_final_${i}" || exit
  mpirun ${LAMMPS_EXEC} ${vars_list}
  if [[ $? -ne 0 ]]; then
    echo "Error running LAMMPS simulation, check the output and log files!"
    exit
  fi
  #         num_files             misorientation                                num_types neigh_cutoff orientation_cutoff a0 structure
  echo "$(find . -name  "*.dump" | wc -l) $(basename "${datafile}" | awk -F'_' '{print $6}' | cut -c 1-5) 2 1.207 0.0 5.454 fcc" > find_grains_input.txt
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
    112)
      echo -e "1 1 1\n1 -1 0\n1 1 -2" >> find_grains_input.txt
      ;;  
    *)  
      echo "ERROR: Unrecognized axis ${axis}"
      exit
      ;;  
  esac
  ls -v *.dump >> find_grains_input.txt
  ${ffg} find_grains_input.txt -e 0 -i 2
  echo "data.txt ${Temp} 0.0 $(awk 'NR==7 {print $2}' "${datafile}") 5.454 fcc" > grain_area_input.txt
  ${cga} grain_area_input.txt --potential-file ${potential_db} -p 1
  tar -cjvf snapshots.tar.bz2 ./*.dump && rm ./*.dump
  )
#done

archive="${datadir}/dir_final_${i}"
echo "Archiving ${archive}"
tar -cjvf "/home/jarinf/archive_tmp/UO2_${axis}_Basak_${mis_dir}_T${Temp}_${radius}_dir${i}_results.tar.bz2" "${archive}" && rm -rf "${archive}" && echo "Successfully archived UO2_${axis}_Basak_${mis_dir}_T${Temp}_${radius}_dir${i}_results!"
