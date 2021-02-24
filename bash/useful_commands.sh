#! /bin/bash

# This command is used when processing grain growth data. To limit to total
# amount of data being saved, this command can be run to get rid of the dump
# files where grain growth has stopped (i.e. the grain has shrunk completely).
for i in 1{0,1}{0,1}/{adp,meam,ternary_eam}/{{20,30,45}degree,sigma7}/T{1050..1400..50}/large_r; do
for i in {5,14,17,21,50}at%/{20,30,45}degree/T{1050..1400..50}/large_r; do
  if [ -d ${i} ]; then
    if [ $(ls "${i}"/dir_*/fitted_values.txt 2> /dev/null | wc -l) -eq 0 ]; then
      echo "Skipping because there are 0 files found for ${i}."
      continue
    fi
    num_to_keep=$(cat ${i}/dir_*/fitted_values.txt | awk -F ':' 'BEGIN{max=0} {if ($2 > max) max = $2} END{print max}' | xargs -I{} qalc -t "ceil({}/10)*10 + 22")
    echo "${i}: ${num_to_keep}"
    for j in dir_{1..5}; do
      (if [ -d "${i}/${j}" ]; then
        cd ${i}/${j}
        # NOTE: double check this gets rid of the files you want with an echo statement! e.g.
        if [[ $(ls -v *0.dump | wc -l 2> /dev/null) -gt 0 ]] && [[ "${num_to_keep}" -lt $(ls *0.dump | wc -l) ]]; then
          echo "Deleting $(qalc -t $(ls *0.dump | wc -l)-${num_to_keep}) files"
          echo "ls -v *0.dump | tail -n +${num_to_keep} | xargs rm"
          # ls -v *0.dump | tail -n +${num_to_keep} | xargs rm
        fi
      fi)
    done
  fi
done

#########################################################################################################################
#########################################################################################################################
#########################################################################################################################
#########################################################################################################################

# This command I used to determine whether to keep the "shortened" list of dump
# files, or the whole list. Displays the total saved space at the end.
tmpfile=$(mktemp 2>/dev/null) || tmpfile=/tmp/input$$
for i in {100,110,111}/{adp,meam,ternary_eam}/{{20,30,45}degree,sigma7}/T{1050..1400..50}/large_r/dir_{1..5}; do
  if [ -d "${i}" ]; then
    (cd "${i}" || continue
    if [ -f "snapshots.7z" ] && [ -f "snapshots_short.7z" ]; then
      if [[ $(stat -c%s "snapshots_short.7z") -lt $(qalc -t -set "exp mode off" "abs(1000000-"$(stat -c%s "snapshots.7z")")") ]]; then
        echo "$(qalc -t -set "exp mode off" "$(stat -c%s "snapshots.7z")"-"$(stat -c%s "snapshots_short.7z")")" >> ${tmpfile}
        echo "snapshots_short.7z is smaller than snapshots.7z by $(numfmt --to iec --format "%.2f" $(qalc -t -set "exp mode off" "$(stat -c%s "snapshots.7z")"-"$(stat -c%s "snapshots_short.7z")")): mv snapshots_short.7z snapshots.7z" # Removing the larger list
        # mv snapshots_short.7z snapshots.7z
      else
        echo "Not enough gains for the shorter list: rm snapshots_short.7z"
        # rm snapshots_short.7z
      fi
    else
      echo "Missing a file in ${i}"
    fi)
  fi
done
echo "Total saved space: $(numfmt --to iec --format "%.2f" $(awk 'BEGIN{sum=0} {sum+=$0} END{print sum}' "${tmpfile}"))"
rm "${tmpfile}"

#########################################################################################################################
#########################################################################################################################
#########################################################################################################################
#########################################################################################################################

# This command checks each directory for a file named "interfaces.7z". If it exists,
# it checks to see if it is corrupted by performing a list command. If so, the
# corrupted archive is removed, and a new archive is created. If the file does
# not exist, the archive is created.
for i in 1{0,1}{0,1}/{adp,meam,ternary_eam}/{{20,30,45}degree,sigma7}/T{1050..1400..50}/large_r/dir_{1..5}; do
  if [ -d ${i} ]; then
    (cd ${i}
    if [ -f "interfaces.7z" ]; then
      7z l interfaces.7z > /dev/null 2>&1
      if [[ "$?" -ne 0 ]]; then
        # rm interfaces.7z
        # 7z a -m0=lzma -mx=9 interfaces.7z *_interface.dat
        echo "rm interfaces.7z"
        echo "7z a -m0=lzma -mx=9 interfaces.7z *_interface.dat"
      fi
    else
      # 7z a -m0=lzma -mx=9 interfaces.7z *_interface.dat
      echo "7z a -m0=lzma -mx=9 interfaces.7z *_interface.dat"
    fi
    )
  fi
done

# This command goes through each directory, and checks to see if the number of
# dump files is equal to the number of files in the snapshots.7z archive. If so,
# a new 7z archive is NOT created. Otherwise, a new one is created, called snapshots_short.7z
for i in {adp,ternary_eam}/{14,17,21}at%/{20,30,45}degree/T{1050..1400..50}/large_r/dir_{1..5}; do
  if [ -d "${i}" ]; then
    (cd ${i}
    if [[ $(ls *0.dump | wc -l) -eq $(7z l -slt snapshots.7z | grep -c 'Path = [^/]*[^.7z]$') ]]; then
      echo "Same number of files in existing archive. Continuing..."
      continue;
    else
      7z a -m0=lzma -mx=9 snapshots_short.7z *0.dump
    fi)
  fi
done

#########################################################################################################################
#########################################################################################################################
#########################################################################################################################
#########################################################################################################################

for angle in {5..180..5}.00; do
  for cutoff in $(seq 0.05 0.05 0.95) .99; do
    echo "LAMMPS_U_110_ternary_eam.dat tmp_${angle}degree_rcut${cutoff} 0.75 ${angle} U bcc 3.542 ${cutoff}" > in.rotate_and_remove
    rotate_and_remove in.rotate_and_remove -o r --no-index
    for T in {900..1400..100}; do
      mpirun lmp_openmpigccfftw -var data_file tmp_${angle}degree_rcut${cutoff}_removed.dat -var T ${T} < equilibrate.in > output_${angle}degree_rcut${cutoff}_T${T}.txt
    done
  done
done
