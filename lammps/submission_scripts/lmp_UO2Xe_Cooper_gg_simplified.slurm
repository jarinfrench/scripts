#! /bin/bash

#SBATCH --nodes=2
#SBATCH --ntasks=256
#SBATCH -t 48:00:00
#SBATCH -p preemptable_q
#SBATCH -A personal
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jarinf@vt.edu
#SBATCH -J UO2Xe_eam_gg_1

#SBATCH --array=28,30-31,33-34,38-40,45,59,62,65-68,70-73,75,106,113,116,118-119,125,128-143%24 # dir_final_1
##SBATCH --array=132-143%24 # dir_final_2
##SBATCH --array=132-143%24 # dir_final_3

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
Xe=(0.25 0.50 0.75 1.00)
T=({3000..2200..400}) #Cooper
DATAFILES=(
${rootdir}/LAMMPS_UO2_100_50x50x5_pCooper_20.00degree_r100A_rcut0.15_annealed.dat
${rootdir}/LAMMPS_UO2_100_50x50x5_pCooper_20.00degree_r75A_rcut0.35_annealed.dat
${rootdir}/LAMMPS_UO2_110_50x35x4_pCooper_20.00degree_r100A_rcut0.05_annealed.dat
${rootdir}/LAMMPS_UO2_110_50x35x4_pCooper_20.00degree_r75A_rcut0.40_annealed.dat
${rootdir}/LAMMPS_UO2_111_35x20x3_pCooper_20.00degree_r100A_rcut0.60_annealed.dat
${rootdir}/LAMMPS_UO2_100_50x50x5_pCooper_45.00degree_r100A_rcut0.35_annealed.dat
${rootdir}/LAMMPS_UO2_100_50x50x5_pCooper_45.00degree_r75A_rcut0.05_annealed.dat
${rootdir}/LAMMPS_UO2_110_50x35x4_pCooper_45.00degree_r100A_rcut0.25_annealed.dat
${rootdir}/LAMMPS_UO2_110_50x35x4_pCooper_45.00degree_r75A_rcut0.30_annealed.dat
${rootdir}/LAMMPS_UO2_111_35x20x3_pCooper_45.00degree_r100A_rcut0.45_annealed.dat
${rootdir}/LAMMPS_UO2_111_35x20x3_pCooper_38.20degree_r100A_rcut0.40_annealed.dat
${rootdir}/LAMMPS_UO2_112_25x30x3_pCooper_45.00degree_r86A_rcut0.55_annealed.dat
) 
p=18
T_idx=$((SLURM_ARRAY_TASK_ID%${#T[@]}))
c_idx=$(($(echo "${SLURM_ARRAY_TASK_ID}/${#T[@]}" | bc) % ${#Xe[@]}))
df_idx=$(($(echo "${SLURM_ARRAY_TASK_ID}/(${#T[@]}*${#Xe[@]})" | bc) % ${#DATAFILES[@]}))
datafile=${DATAFILES[${df_idx}]}
c=${Xe[${c_idx}]}
axis=$(basename "${datafile}" | awk -F'_' '{print $3}')
radius=$(basename "${datafile}" | awk -F'_' '{print $7}')
if [[ ${radius} == "r100A" ]] || [[ ${radius} == "r86A" ]]; then
  radius_dir="large_r"
elif [[ ${radius} == "r75A" ]]; then
  radius_dir="medium_r"
else
  echo "Unrecognized grain radius: ${radius} (basename = ${datafile})"
  exit
fi
misorientation=$(basename "${datafile}" | awk -F'_' '{print $6}')
if [[ ${misorientation} == "38.20degree" ]]; then
  mis_dir="sigma7"
else
  mis_dir=$(echo "${misorientation}" | cut -c 3-5 --complement) # cut out the .00 from the file name, this gives us the misorientation directory (mis_dir)
fi

if [[ ${c} == "0.25" ]]; then
  if [[ ${axis} == "100" ]]; then
    if [[ ${misorientation} == "20.00degree" ]]; then
      if [[ ${radius} == "r100A" ]]; then
        timesteps=(1000000 1400000 4000000)
      else # ${radius} == "r75A" 
        timesteps=(500000 800000 2000000)
      fi
    elif [[ ${misorientation} == "45.00degree" ]]; then
      if [[ ${radius} == "r100A" ]]; then
        timesteps=(1000000 1400000 4000000)
      else # ${radius} == "r75A" 
        timesteps=(1000000 1400000 2000000)
      fi
    else # ${misorientation} == "38.20degree"
      echo "wrong axis"
      exit
    fi
  elif [[ ${axis} == "110" ]]; then
    if [[ ${misorientation} == "20.00degree" ]]; then
      if [[ ${radius} == "r100A" ]]; then
        timesteps=(1200000 1900000 4000000)
      else # ${radius} == "r75A" 
        timesteps=(1200000 1900000 4000000)
      fi
    elif [[ ${misorientation} == "45.00degree" ]]; then
      if [[ ${radius} == "r100A" ]]; then
        timesteps=(1500000 2500000 4000000)
      else # ${radius} == "r75A" 
        timesteps=(1500000 2500000 4000000)
      fi
    else # ${misorientation} == "38.20degree"
      echo "wrong axis"
      exit
    fi
  elif [[ ${axis} == "111" ]]; then
    if [[ ${misorientation} == "20.00degree" ]]; then
      if [[ ${radius} == "r100A" ]]; then
        timesteps=(1200000 1400000 2500000)
      else # ${radius} == "r75A" 
        echo "wrong radius"
        exit
      fi
    elif [[ ${misorientation} == "45.00degree" ]]; then
      if [[ ${radius} == "r100A" ]]; then
        timesteps=(1200000 1800000 3000000)
      else # ${radius} == "r75A" 
        echo "wrong radius"
        exit
      fi
    else # ${misorientation} == "38.20degree"
       timesteps=(1200000 1800000 3000000)
    fi
  else # ${axis} == "112"
    if [[ ! ${misorientation} == "45.00degree" ]] && [[ ! ${radius} == "r86A" ]]; then
      echo "misorientation or radius wrong"
      exit
    else
      timesteps=(1500000 2200000 3000000)
    fi
  fi
elif [[ ${c} == "0.50" ]]; then
  if [[ ${axis} == "100" ]]; then
    if [[ ${misorientation} == "20.00degree" ]]; then
      if [[ ${radius} == "r100A" ]]; then
        timesteps=(1800000 2200000 4000000)
      else # ${radius} == "r75A" 
        timesteps=(1300000 1600000 3000000)
      fi
    elif [[ ${misorientation} == "45.00degree" ]]; then
      if [[ ${radius} == "r100A" ]]; then
        timesteps=(1800000 2200000 4000000)
      else # ${radius} == "r75A" 
        timesteps=(1800000 2200000 3000000)
      fi
    else # ${misorientation} == "38.20degree"
      echo "wrong axis"
      exit
    fi
  elif [[ ${axis} == "110" ]]; then
    if [[ ${misorientation} == "20.00degree" ]]; then
      if [[ ${radius} == "r100A" ]]; then
        timesteps=(2000000 3000000 4000000)
      else # ${radius} == "r75A" 
        timesteps=(2000000 3000000 4000000)
      fi
    elif [[ ${misorientation} == "45.00degree" ]]; then
      if [[ ${radius} == "r100A" ]]; then
        timesteps=(2300000 3300000 4000000)
      else # ${radius} == "r75A" 
        timesteps=(2300000 3300000 4000000)
      fi
    else # ${misorientation} == "38.20degree"
      echo "wrong axis"
      exit
    fi
  elif [[ ${axis} == "111" ]]; then
    if [[ ${misorientation} == "20.00degree" ]]; then
      if [[ ${radius} == "r100A" ]]; then
        timesteps=(2000000 2200000 3300000)
      else # ${radius} == "r75A" 
        echo "wrong radius"
        exit
      fi
    elif [[ ${misorientation} == "45.00degree" ]]; then
      if [[ ${radius} == "r100A" ]]; then
        timesteps=(2000000 2600000 4000000)
      else # ${radius} == "r75A" 
        echo "wrong radius"
        exit
      fi
    else # ${misorientation} == "38.20degree"
       timesteps=(2000000 2600000 4000000)
    fi
  else # ${axis} == "112"
    if [[ ! ${misorientation} == "45.00degree" ]] && [[ ! ${radius} == "r86A" ]]; then
      echo "misorientation or radius wrong"
      exit
    else
      timesteps=(2300000 3000000 4000000)
    fi
  fi
elif [[ ${c} == "0.75" ]]; then
  if [[ ${axis} == "100" ]]; then
    if [[ ${misorientation} == "20.00degree" ]]; then
      if [[ ${radius} == "r100A" ]]; then
        timesteps=(2600000 3000000 4000000)
      else # ${radius} == "r75A" 
        timesteps=(2100000 2400000 4000000)
      fi
    elif [[ ${misorientation} == "45.00degree" ]]; then
      if [[ ${radius} == "r100A" ]]; then
        timesteps=(2600000 3000000 4000000)
      else # ${radius} == "r75A" 
        timesteps=(2600000 3000000 4000000)
      fi
    else # ${misorientation} == "38.20degree"
      echo "wrong axis"
      exit
    fi
  elif [[ ${axis} == "110" ]]; then
    if [[ ${misorientation} == "20.00degree" ]]; then
      if [[ ${radius} == "r100A" ]]; then
        timesteps=(3000000 4000000 4000000)
      else # ${radius} == "r75A" 
        timesteps=(3000000 4000000 4000000)
      fi
    elif [[ ${misorientation} == "45.00degree" ]]; then
      if [[ ${radius} == "r100A" ]]; then
        timesteps=(3100000 4000000 4000000)
      else # ${radius} == "r75A" 
        timesteps=(3100000 4000000 4000000)
      fi
    else # ${misorientation} == "38.20degree"
      echo "wrong axis"
      exit
    fi
  elif [[ ${axis} == "111" ]]; then
    if [[ ${misorientation} == "20.00degree" ]]; then
      if [[ ${radius} == "r100A" ]]; then
        timesteps=(2800000 3000000 4000000)
      else # ${radius} == "r75A" 
        echo "wrong radius"
        exit
      fi
    elif [[ ${misorientation} == "45.00degree" ]]; then
      if [[ ${radius} == "r100A" ]]; then
        timesteps=(2800000 3400000 4000000)
      else # ${radius} == "r75A" 
        echo "wrong radius"
        exit
      fi
    else # ${misorientation} == "38.20degree"
       timesteps=(2800000 3400000 4000000)
    fi
  else # ${axis} == "112"
    if [[ ! ${misorientation} == "45.00degree" ]] && [[ ! ${radius} == "r86A" ]]; then
      echo "misorientation or radius wrong"
      exit
    else
      timesteps=(3100000 3800000 4000000)
    fi
  fi
else # ${c} == "1.00"
  if [[ ${axis} == "100" ]]; then
    if [[ ${misorientation} == "20.00degree" ]]; then
      if [[ ${radius} == "r100A" ]]; then
        timesteps=(3400000 4000000 4000000)
      else # ${radius} == "r75A" 
        timesteps=(3000000 4000000 4000000)
      fi
    elif [[ ${misorientation} == "45.00degree" ]]; then
      if [[ ${radius} == "r100A" ]]; then
        timesteps=(3400000 4000000 4000000)
      else # ${radius} == "r75A" 
        timesteps=(3400000 4000000 4000000)
      fi
    else # ${misorientation} == "38.20degree"
      echo "wrong axis"
      exit
    fi
  elif [[ ${axis} == "110" ]]; then
    if [[ ${misorientation} == "20.00degree" ]]; then
      if [[ ${radius} == "r100A" ]]; then
        timesteps=(4000000 4000000 4000000)
      else # ${radius} == "r75A" 
        timesteps=(4000000 4000000 4000000)
      fi
    elif [[ ${misorientation} == "45.00degree" ]]; then
      if [[ ${radius} == "r100A" ]]; then
        timesteps=(4000000 4000000 4000000)
      else # ${radius} == "r75A" 
        timesteps=(4000000 4000000 4000000)
      fi
    else # ${misorientation} == "38.20degree"
      echo "wrong axis"
      exit
    fi
  elif [[ ${axis} == "111" ]]; then
    if [[ ${misorientation} == "20.00degree" ]]; then
      if [[ ${radius} == "r100A" ]]; then
        timesteps=(3600000 4000000 4000000)
      else # ${radius} == "r75A" 
        echo "wrong radius"
        exit
      fi
    elif [[ ${misorientation} == "45.00degree" ]]; then
      if [[ ${radius} == "r100A" ]]; then
        timesteps=(3600000 4000000 4000000)
      else # ${radius} == "r75A" 
        echo "wrong radius"
        exit
      fi
    else # ${misorientation} == "38.20degree"
       timesteps=(3600000 4000000 4000000)
    fi
  else # ${axis} == "112"
    if [[ ! ${misorientation} == "45.00degree" ]] || [[ ! ${radius} == "r86A" ]]; then
      echo "misorientation or radius wrong"
      exit
    else
      timesteps=(4000000 4000000 4000000)
    fi
  fi
fi

dt=${timesteps[${T_idx}]}
Temp=${T[${T_idx}]}

datadir=${rootdir}/${c}%Xe/${axis}/${mis_dir}/${potential}/T${Temp}/${radius_dir}
echo "${datadir} being processed" >> /home/jarinf/projects/uo2/impurity_effect/grain_growth/slurm_log.txt 

i=1
(
mkdir -p "${datadir}/dir_final_${i}"
SEED=${RANDOM}
impurityfile=$(basename "${datafile}" .dat | awk -v Xe="${c}" -F'UO2' '{print $1"UO2+"Xe"Xe"$2".dat"}')
echo "${datafile} ${datadir}/${impurityfile} UO2 $(printf '%.4f\n' "$(echo "scale=4; ${c}/100" | bc)") 1 ${RANDOM} 5.454" | tee "input.tmp${SLURM_ARRAY_TASK_ID}" >> "${rootdir}/generate_impurities_input.log"
${impgen} "input.tmp${SLURM_ARRAY_TASK_ID}" -d s
vars_list="-var SEED ${SEED} -var file ${datadir}/${impurityfile} -var T ${Temp} -var timesteps ${dt} -var cXe ${c} -var newfile $(echo "$(basename ${impurityfile} .dat)_annealed.dump")"
cd "${datadir}/dir_final_${i}" || exit
mpirun ${LAMMPS_EXEC} ${vars_list} < ${INFILE} > "output_${Temp}.txt"
#         num_files             misorientation                                num_types neigh_cutoff orientation_cutoff a0 structure
echo "$(ls *.dump | wc -l) $(basename "${datafile}" | awk -F'_' '{print $6}' | cut -c 1-5) 3 1.207 0.0 5.454 fcc" > find_grains_input.txt
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
ls -v ./*.dump >> find_grains_input.txt
${ffg} find_grains_input.txt -e 0 -i 2
echo "data.txt ${Temp} $(printf '%.3f\n' "$(echo "scale=4; ${c}/100" | bc)") $(awk 'NR==7 {print $2}' "${datafile}") 5.454 fcc" > grain_area_input.txt
${cga} grain_area_input.txt --potential-file ${potential_db} -p ${p}
tar -cjvf snapshots.tar.bz2 ./*.dump && rm ./*.dump
)

archive="${datadir}/dir_final_${i}"
echo "Archiving ${archive}"
tar -cjvf "/home/jarinf/archive_tmp/UO2+${c}Xe_${axis}_${mis_dir}_T${Temp}_${radius}_dir${i}_results.tar.bz2" "${archive}" && rm -rf "${archive}" && echo "Successfully archived UO2+${c}Xe_${axis}_${mis_dir}_T${Temp}_${radius}_dir${i}_results!"

