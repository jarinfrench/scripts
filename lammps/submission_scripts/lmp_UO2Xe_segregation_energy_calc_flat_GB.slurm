#! /bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH -t 36:00:00
#SBATCH -p normal_q
#SBATCH -A GB_Fuels
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jarinf@vt.edu
#SBATCH -J UO2_XePu_seg

module reset
module load foss/2020a

LAMMPS_EXEC=/home/jarinf/LAMMPS/lammps/tinkercliffs_mpi_build_foss/lmp_29Oct2020
impgen=/home/jarinf/scripts/tinkercliffs_foss_2020a_executables/generate_impurities
DATAFILES=($(ls *_wrapped.lmp))

for file in ${DATAFILES[@]}; do
  axis=$(echo ${file} | awk -F'_' '{print $2}')
  csl=$(echo ${file} | awk -F'_' '{print $3}')
  while read -u 9 id; do
    echo "${file} tmp_${id}.dat UO2 ${id} 5.454" | tee tmp_impgen.txt >> impurity_generation_inputs.log
    ${impgen} tmp_impgen.txt -d s --no-molecule-removal
  
    # Run the Basak case
    vars="-var id ${id} -var potential Basak -var axis ${axis} -var mis ${csl} -var datafile tmp_${id}.dat -in segregation_Xe.in -sc output_Basak+Xe_${axis}_${csl}_${id}.txt"
    mpirun ${LAMMPS_EXEC} ${vars}
  
    # Run the Cooper case for Xe
    vars="-var id ${id} -var potential Cooper -var axis ${axis} -var mis ${csl} -var datafile tmp_${id}.dat -in segregation_Xe.in -sc output_Cooper+Xe_${axis}_${csl}_${id}.txt"
    mpirun ${LAMMPS_EXEC} ${vars}
  
    # Run the Cooper case for Pu
    vars="-var id ${id} -var potential Cooper -var axis ${axis} -var mis ${csl} -var datafile tmp_${id}.dat -in segregation_Pu.in -sc output_Cooper+Pu_${axis}_${csl}_${id}.txt"
    mpirun ${LAMMPS_EXEC} ${vars}
  
    rm tmp_${id}.dat
  done 9< ${file%.lmp}_atom_ids.txt
done

tar -cjvf outputs_flat_GB_U_replace.tar.bz2 output_*.txt && rm output_*.txt
