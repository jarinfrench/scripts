#! /bin/bash

# Given an input script, lmp.pbs for example, variables can be passed to
# qsub using the -v command, i.e. qsub -v T=10 lmp.pbs.  The variable T in lmp.pbs
# would then be replaced with the value 10 everywhere that $T appears.


NUM_CORES=16 # Determined after doing a scaling study

read -p "Please enter the simulation temperature: " TEMPERATURE
read -p "Please enter the initial position data file: " INITIAL_POSITION
read -p "How many timesteps will be run: " NUM_TIMESTEPS
read -p "Please enter the element: " ELEMENT

ELEMENT=$(echo "${ELEMENT}" | tr '[:upper:]' '[:lower:]')
if [[ "${ELEMENT}" =~ uo2 ]]; then
  ELEMENT=$(echo "${ELEMENT}" | tr '[:lower:]' '[:upper:]')
else
  ELEMENT="$(tr '[:lower:]' '[:upper:]' <<< ${ELEMENT:0:1})${ELEMENT:1}"
fi

LAMMPS_INPUT="${ELEMENT}_MSD_gg.in"

ARGS="-var SEED \`bash -c 'echo $``RANDOM'\` -var T ${TEMPERATURE} -var data_file ${INITIAL_POSITION} -var timesteps $NUM_TIMESTEPS"
read -p "Please specify either ARC or INL: " HPC
HPC=`echo ${HPC} | awk '{print tolower($0)}'`
if [[ !("${HPC}" =~ ^(inl|arc)$) ]]; then
  echo "${HPC} not recognized.  Please enter either INL or ARC."
  exit 4
fi

# Note that the multiplier is specific to UO2 systems using the Basak potential.  Different systems will have different multipliers.
# This could probably be updated by making the multiplier be the walltime per timestep per atom, and prompt for the number of atoms in the simulation. TODO!
if [[ "${HPC}" == "inl" ]]; then
  NUM_PROCS=24 # Trying to run HPC with all 36 processors results in segfaults and BUS errors
  NUM_CPUS=36
  MULTIPLIER=7.86e-6 # walltime per timestep for UO2 in Falcon
elif [[ "${HPC}" == "arc" ]]; then
  NUM_PROCS=32
  NUM_CPUS=32
  MULTIPLIER=7.356e-6 # walltime per timestep for UO2 in Cascades
else
  echo "Something went wrong.  HPC system not recognized."
  exit 5
fi

# We add the extra two just in case something weird happens during the simulation that takes longer than expected.
WALLTIME=`python -c "from math import ceil; print '%.0f'%(ceil(${NUM_TIMESTEPS}*${MULTIPLIER})+2)"`
echo "Estimated simulation time is ${WALLTIME} hours"

# Now create the script
if [[ "${HPC}" == "inl" ]]; then # if INL
  echo "#! /bin/bash

#PBS -j oe
#PBS -l select=${NUM_CORES}:ncpus=${NUM_CPUS}:mpiprocs=${NUM_PROCS}
#PBS -l walltime=${WALLTIME}:00:00
#PBS -P neams

#PBS -M jarinf@vt.edu
#PBS -m bea
#PBS -N T${TEMPERATURE}

source /etc/profile.d/modules.sh
module load OpenMPI

cd $``PBS_O_WORKDIR
mkdir dir_$``n
cp lmp.pbs base_vals.txt ${INITIAL_POSITION} ${LAMMPS_INPUT} dir_$``n
cd dir_$``n

mpirun /home/frenjari/projects/lammps/lammps-17Nov16/src/lmp_openmpigccfftw_falcon ${ARGS} < ${LAMMPS_INPUT} > output.txt

cd ..
exit;" > lmp.pbs
elif [[ "${HPC}" == "arc" ]]; then
  echo "#! /bin/bash

#PBS -j oe
#PBS -l nodes=${NUM_CORES}:ppn=${NUM_PROCS}
#PBS -l walltime=${WALLTIME}:00:00
#PBS -q normal_q
#PBS -A GB_Fuels
#PBS -W group_list=cascades

#PBS -M jarinf@vt.edu
#PBS -m bea
#PBS -N T${TEMPERATURE}

module purge
module load gcc openmpi

cd $``PBS_O_WORKDIR
mkdir dir_$``n
cp lmp.pbs base_vals.txt ${INITIAL_POSITION} ${LAMMPS_INPUT} dir_$``n
cd dir_$``n

mpirun -np $``PBS_NP /home/jarinf/LAMMPS/lammps-17Nov16/src/lmp_openmpigccfftw_cascades ${ARGS} < ${LAMMPS_INPUT} > output.txt

cd ..
exit;" > lmp.pbs
fi
