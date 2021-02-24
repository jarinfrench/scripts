#! /bin/bash

while [[ ${scheduler,,} != "pbs" ]] && [[ ${scheduler,,} != "slurm" ]] || [[ -z ${scheduler} ]]; do
  read -p "Enter the scheduler (pbs|slurm): " scheduler
done

if [[ ${scheduler,,} == "pbs" ]]; then
  echo "Not implemented"
  exit 0
fi

while [[ -z ${nodes} ]]; do
  read -p "Enter the number of nodes: " nodes
done

while [[ -z ${tasks} ]]; do
  read -p "Enter the number of processors per node: " tasks
done

while [[ -z ${walltime} ]] || [[ ! ${walltime} =~ [0-9][0-9]:[0-5][0-9]:[0-5][0-9] ]]; do
  read -p "Enter the walltime (hr:mi:se): " walltime
done

while [[ -z ${job_name} ]]; do
  read -p "Enter job name: " job_name
done

if [[ ${scheduler,,} == "slurm" ]]; then
  while [[ -z ${account} ]]; do # TODO: also add range checking...
    read -p '''Select account:
1 - gb_fuels
2 - fecr_bai
3 - conc_alloys
4 - oxides_1
5 - ris_fracture

Account: ''' account
  done

  directive="#SBATCH "
  dest_file="${job_name}.${scheduler,,}"

  echo -e "#! /bin/bash\n" >> ${dest_file}
  echo "${directive} --nodes=${nodes}" >> ${dest_file}
  echo "${directive} --ntasks=$((${tasks}*${nodes}))" >> ${dest_file}
  echo "${directive} -t ${walltime}" >> ${dest_file}
  echo "${directive} -p normal_q" >> ${dest_file}
  echo "${directive} -A ${account}" >> ${dest_file}
  echo "${directive} --mail-type=END,FAIL" >> ${dest_file}
  echo "${directive} --mail-user=jarinf@vt.edu" >> ${dest_file}
  echo "${directive} -J ${job_name}" >> ${dest_file}
  echo -e "\nmodule purge" >> ${dest_file}
  echo "module load <specify modules here>" >> ${dest_file}
fi
