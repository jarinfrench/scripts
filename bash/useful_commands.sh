#! /bin/bash

# This command is used when processing grain growth data. To limit to total
# amount of data being saved, this command can be run to get rid of the dump
# files where grain growth has stopped (i.e. the grain has shrunk completely).
#for i in {5,14,17,21,50}at%/{20,30,45}degree/T{1050..1400..50}/large_r; do
for i in 1{0,1}{0,1}/{adp,{m,ternary_}eam}/{{20,30,45}degree,sigma7}/T{900..1400..50}/large_r; do
  if [ -d "${i}" ]; then # directory exists
    if [ $(ls "${i}"/dir_*/slope_calc.txt 2> /dev/null | wc -l) -eq 0 ]; then # counts the number of slope_calc.txt files found in subdirs of ${i}
      slope_calc_found=0
    else
      slope_calc_found=1
    fi
    if [ $(ls "${i}"/dir_*/fitted_values.txt 2> /dev/null | wc -l) -eq 0 ]; then # similarly for fitted_values.txt files
      fitted_values_found=0
    else
      fitted_values_found=1
    fi
    if [ ${slope_calc_found} -eq 0 ] && [ ${fitted_values_found} -eq 0 ]; then
      echo "No slope_calc.txt or fitted_values.txt found in ${i}" >> tmp.none_found.txt
      continue
    fi
    # the ceil({}/10)*10 + 12 comes from the fact we want a factor of 10 files kept. The +12 keeps 12 values past the calculated factor of 10
    if [ ${slope_calc_found} -eq 1 ]; then # prefer slope_calc.txt files
      # the first if condition in the awk statement takes care of slope_calc.txt files where the whole line matches the first field ($0 == $1). If this is the case, keep ALL dump files
      num_to_keep=$(cat ${i}/dir_*/slope_calc.txt | awk -F ':' 'BEGIN{max=0} {if ($0 == $1) max = 10000; else if ($2 > max) max = $2} END{print max}' | xargs -I{} qalc -t "ceil({}/10)*10 + 12")
    elif [ ${fitted_values_found} -eq 1 ]; then
      num_to_keep=$(cat ${i}/dir_*/fitted_values.txt | awk -F ':' 'BEGIN{max=0} {if ($2 > max) max = $2} END{print max}' | xargs -I{} qalc -t "ceil({}/10)*10 + 12")
    fi
    # echo "${i}: ${num_to_keep}" >> tmp.kept_values.txt
    for j in dir_{1..10}; do
      (if [ -d "${i}/${j}" ]; then
        cd ${i}/${j}
        if [[ $(ls *0.dump 2> /dev/null | wc -l ) -gt 12 ]] && [[ "${num_to_keep}" -lt $(ls *0.dump 2> /dev/null | wc -l ) ]]; then
          echo "Will delete $(qalc -t $(ls *0.dump | wc -l)-${num_to_keep}) out of $(ls *0.dump | wc -l) files"
          echo "ls -v *0.dump | tail -n +${num_to_keep} | xargs rm"
          # ls -v *0.dump | tail -n +${num_to_keep} | xargs rm
        else
          echo "$(ls *0.dump 2> /dev/null| wc -l || echo 0) dump files found (keeping ${num_to_keep}) in ${i}/${j}"
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

#########################################################################################################################
#########################################################################################################################
#########################################################################################################################
#########################################################################################################################
# This command travels through the U grain growth directories, and extracts
# the first and last dump file from the snapshots.7z archive if it exists.
# if it doesn't exist, it looks for the uncompressed files. If found, the
# script then runs the command track_atoms on those files to generate a
# file that can be used by tecplot to examine the rotational behavior of the
# simulations

for i in {100,110,111}/{adp,meam,ternary_eam}/{{20,30,45}degree,sigma7}/T{1050..1400..50}/large_r/dir_{1..5}; do
  if [[ -d "${i}" ]]; then
    (cd ${i}
    if [[ -f "snapshots.7z" ]]; then
      last_file=$(7z l snapshots.7z | awk '/------------------- ----- ------------ ------------  ------------------------/{n+=1; r+=1} n % 2 == 1 && ! /------------------- ----- ------------ ------------  ------------------------/  {if (r == 1) {print $6; r -= 1} else {print $5}}' | sort -V | tail -n1)
      7z x -mmt=off ./snapshots.7z -aos "${last_file}" > /dev/null
      7z x -mmt=off ./snapshots.7z -aos "0.dump" > /dev/null
    else
      echo -e "\033[0;31msnapshots.7z not found in ${i}\033[0m"
      if test -n "$(shopt -s nullglob; echo *.dump)"; then
        last_file=$(ls -v *.dump | tail -n1)
      else
        echo -e "\033[0;31mNo dump files found in ${i}\033[0m"
        break
      fi
    fi
    track_atoms 0.dump "${last_file}"
    )
  fi
done

#########################################################################################################################
#########################################################################################################################
#########################################################################################################################
#########################################################################################################################
# This command find all the dir directories and does the following:
#   If a snapshots.7z archive exists, it is extracted
#   If not, it is created from the existing *.dump files
#   If an interfaces.7z archive exists, it is extracted
#   If not, it is created from the existing *_interface.dat files
#   Creates a file for the tracked atoms in the central strip
#   Creates two files for the tracked atoms in a radius centered on the grain
#   If a snapshots.7z archive exists, removes all the dump files
#   Else creates the archive, then removes the files
#   If an interfaces.7z archive exists, removes all the *_interface.dat files
#   Else creates the archive, then removes the files

for i in $(fd dir -E plots -t d); do # for each dir found that is a directory, excluding the plots directory
  (
  if [[ -d ${i} ]]; then # if the found value is a directory (double check!)
    echo "Entering ${i}"
    cd ${i}; # change to the directory
    # if [[ -f snapshots.7z ]]; then # if a snapshots.7z file exists
    #   extract snapshots.7z; # extract it
    # else
    #   7za snapshots.7z *.dump; # otherwise create the snapshots.7z file
    # fi;
    # if [[ -f interfaces.7z ]]; then # if an interfaces.7z file exists
    #   extract interfaces.7z; # extract it
    # else
    #   7za interfaces.7z *_interface.dat; # otherwise, make it
    # fi;
    if [[ -f "final_strip_tracked.dat" && -f "0_tracked.dat" ]]; then
      continue;
    fi
    extract snapshots.7z interfaces.7z

    track_atoms 0.dump $(ls -v *.dump | tail -n 1) -o final_strip_tracked.dat; # track the atoms via a strip

    bounds=$(awk -F ' ' '{if ($1 == "orthogonal") {print $4,$8,$5,$9,$6,$10; exit}}' output.txt | tr -d '()');
    bounds_tmp=(${bounds});
    if [[ ${#bounds_tmp[@]} -ne 6 ]]; then
      echo "Bounds variable not set correctly: bounds = ${bounds}";
    fi;

    dumps=($(ls -v *.dump));
    interfaces=($(ls -v *_interface.dat));

    if [[ -f slope_calc.txt ]]; then
      idx=$(tail -n1 slope_calc.txt | awk -F':' '{print $2}');
    elif [[ -f fitted_values.txt ]]; then
      idx=$(tail -n1 fitted_values.txt | awk -F':' '{print $2}');
    else
      echo "Unable to determine grain growth index";
      continue;
    fi

    if [[ ${idx} -lt 50 ]]; then
      idx=$((${idx}-10))
    elif [[ ${idx} -lt 100 ]]; then
      idx=$((${idx}-20))
    elif [[ ${idx} -lt 150 ]]; then
      idx=$((${idx}-30))
    elif [[ ${idx} -lt 200 ]]; then
      idx=$((${idx}-40))
    else
      idx=$((${idx}-50))
    fi

    dump_num=$(echo "${dumps[${idx}]%.*}")
    idx2=0
    for j in "${interfaces[@]%_*}"; do
      if [[ $((${dump_num}-${j})) -lt 0 ]]; then
        idx2=$((${idx2}-1));
        break;
      else
        idx2=$((${idx2}+1));
      fi;
    done;
    final_center=($(grain_center_calculator.py ${interfaces[${idx2}]} -B ${bounds} -g 2 -pt));
    track_atoms 0.dump ${dumps[${idx}]} --r-type c --rhi 30 --cx ${final_center[0]} --cy ${final_center[1]} --cz ${final_center[2]};
    rm *.dump *_interface.dat
    # if [[ -f snapshots.7z ]]; then
    #   rm *.dump;
    # else
    #   7za snapshots.7z *.dump;
    #   rm *.dump;
    # fi;
    # if [[ -f interfaces.7z ]]; then
    #   rm *_interface.dat;
    # else
    #   7za interfaces.7z *_interface.dat;
    #   rm *_interface.dat;
    # fi;
  fi;
  );
done

# This command creates the slope_calc.txt file (via mobility_plots.py) if it does
# not exist, and if it does exist, it makes sure the file ends in a new line
for i in {100,110,111}/{adp,meam,ternary_eam}/{{20,30,45}degree,sigma7}/T{900..1400..50}/large_r/dir_{1..10}; do
  if [[ -d "${i}" ]]; then
    if [[ ! -f "${i}/slope_calc.txt" ]]; then
      (cd ${i};
      mobility_plots.py area_data.txt -F
      );
    elif [[ -s "${i}/slope_calc.txt" && -z "$(tail -c 1 "${i}/slope_calc.txt")" ]]; then
      (
      cd ${i};
      sed -i -e :a -e '/^\n*$/{$d;N;ba' -e '}' slope_calc.txt
      );
    fi;
  fi;
done

cd '/media/jarinf/Research Backup/Research1_backup/U/grain_growth/gamma/'
for i in $(fd dir -E plots -t d); do # for each dir found that is a directory, excluding the plots directory
  (
  basename=$(echo "${i}" | awk -F'/' '{print $2"_"$1"_"$3"_"$4"_"$6}')
  if [[ -f "${i}/snapshots.7z" ]]; then
    if [[ -f "${i}/slope_calc.txt" ]]; then
      idx=$(tail -n1 "${i}/slope_calc.txt" | awk -F':' '{print $2}');
    else
      echo "Unable to determine grain growth index in ${i}";
      continue;
    fi
    if [[ ${idx} -lt 50 ]]; then
      idx=$((${idx}-10))
    elif [[ ${idx} -lt 100 ]]; then
      idx=$((${idx}-20))
    elif [[ ${idx} -lt 150 ]]; then
      idx=$((${idx}-30))
    elif [[ ${idx} -lt 200 ]]; then
      idx=$((${idx}-40))
    else
      idx=$((${idx}-50))
    fi

    dumps=($(7z l "${i}/snapshots.7z" | grep 0.dump | awk '{print $NF}' | sort -g)) # get the sorted list of dump files
    if [[ "${idx}" -gt "${#dumps[@]}" ]]; then
      echo "Specified dump file does not exist in ${i}"
      continue
    fi
    7z x -mmt=off -aos "${i}/snapshots.7z" -o/media/jarinf/Research2/working 0.dump "${dumps[${idx}]}" # extract the specific files we need
  else
    echo "snapshots.7z not found in ${i}"
    continue
  fi

  if [[ -f "${i}/output.txt" ]]; then
    bounds=$(awk -F ' ' '{if ($1 == "orthogonal") {print $4,$8,$5,$9,$6,$10; exit}}' "${i}/output.txt" | tr -d '()'); # get the bounds from the LAMMPS output file
    bounds_tmp=(${bounds}); # make sure there are 6 elements!
    if [[ ${#bounds_tmp[@]} -ne 6 ]]; then
    echo "Bounds variable not set correctly: bounds = ${bounds}";
    fi;
  else
    echo "output.txt file not found in ${i}. Unable to determine system bounds"
    continue;
  fi

  if [[ -f "${i}/interfaces.7z" ]]; then
    interfaces=($(7z l "${i}/interfaces.7z" | grep interface.dat | awk '{print $NF}' | sort -g)) # get the sorted list of interface files
    dump_num=$(echo "${dumps[${idx}]%.*}") # get the number of the dump file (e.g. 100000 from 100000.dump)
    idx2=0
    for j in "${interfaces[@]%_*}"; do
      if [[ $((${dump_num}-${j})) -lt 0 ]]; then
        idx2=$((${idx2}-1));
        break;
      else
        idx2=$((${idx2}+1));
      fi;
    done;

    # extract the specific interface file we need
    7z x -mmt=off -aos "${i}/interfaces.7z" -o/media/jarinf/Research2/working "${interfaces[${idx2}]}"
    final_center=($(grain_center_calculator.py "/media/jarinf/Research2/working/${interfaces[${idx2}]}" -B ${bounds} -g 2 -pt));
    track_atoms /media/jarinf/Research2/working/0.dump --r-type c --rhi 30 --cx ${final_center[0]} --cy ${final_center[1]} --cz ${final_center[2]} --print-atom-ids > /media/jarinf/Research2/working/atom_ids.txt
  else
    echo "interfaces.7z file not found in ${i}"
    continue
  fi
  calculate_displacement /media/jarinf/Research2/working/0.dump "/media/jarinf/Research2/working/${dumps[${idx}]}" -o "/media/jarinf/Research2/working/${basename}_displacement_data_circle.dat" -f /media/jarinf/Research2/working/atom_ids.txt
  rm /media/jarinf/Research2/working/*.dump /media/jarinf/Research2/working/*_interface.dat /media/jarinf/Research2/working/atom_ids.txt
  cp "/media/jarinf/Research2/working/${basename}_displacement_data.dat" "${i}/0to${dumps[${idx}]%.*}_displacement_data"
  );
done

# This is currently set up for basak 112 45degree
# for i in {2000..2400..50} {2500..3300..100} 3050 3150 ; do
# for i in {2100..2400..50} {2500..3000..100} 3050; do # EAM 110 45 degree
for i in {2000..3300..100}; do
  for j in {large,medium,small,tiny}_r{,/big}/dir_{1..12}; do
    if [[ -d "T${i}/${j}" ]]; then
      cd "T${i}/${j}"
    else
      continue
    fi
    7za old_results.7z *input* area_data.txt data.txt fitted_values.txt *base_vals*.txt slope_calc.txt
    rm *input* area_data.txt data.txt fitted_values.txt *base_vals*.txt slope_calc.txt
    if [[ -f interfaces.7z ]]; then
      rm interfaces.7z
    fi
    if [[ -d interfaces ]]; then
      rm -rf interfaces
    fi
    if [[ -f snapshots.7z ]]; then
      extract snapshots.7z
    fi
    echo "$(ls *.dump | wc -l) 45.00 2 1.207 0.0 5.453 fcc" > find_grains_input.txt
    echo "1 1 1" >> find_grains_input.txt
    echo "1 -1 0" >> find_grains_input.txt
    echo "1 1 -2" >> find_grains_input.txt
    ls -v *.dump >> find_grains_input.txt
    find_grains find_grains_input.txt -e 10 -i 2
    echo "data.txt ${i} 0.0 40.071200 5.453 fcc" > grain_area_input.txt
    calculate_grain_area grain_area_input.txt -p 1
    rm *.dump
    cd "${OLDPWD}"
  done
done

for i in {2000..3300..100}; do
  for j in {large,medium,small,tiny}_r{,/big}/dir_{1..12}; do
    if [[ -d "T${i}/${j}" ]]; then
      cd "T${i}/${j}"
    else
      continue
    fi
    if [[ -f snapshots.7z ]]; then
      extract snapshots.7z
    fi
    find_grains find_grains_input.txt -e 10 -i 2
    calculate_grain_area grain_area_input.txt -p 1
    rm *.dump
    cd "${OLDPWD}"
  done
done

# written specifically to calculate the vacancy diffusion of various U(Mo|Xe) systems using the EAM potential and the ADP potentials
# See ${RESEARCH_HOME}/U/vacancy_diffusion
for i in {0.05,0.10,0.14,0.17,0.21,0.30,0.5}; do # remember to include 0 if doing both Xe and Mo
  # for k in {0,0.001,0.005,0.01,0.02,0.05}; do
    for j in {1050,1200,1400}; do
      # if [[ $(qalc -t ${k} == 0) -eq 1 ]] && [[ $(qalc -t ${i} != 0) -eq 1 ]]; then
        nU=$(awk -v co=${i} -v T=${j} '{if ($1 == co && $2 == T) print $3}' numbers.txt); # these two are for the ADP potential
        nMo=$(awk -v co=${i} -v T=${j} '{if ($1 == co && $2 == T) print $4}' numbers.txt);
        # nU=$(awk -v co=${i} -v T=${j} '{if ($1 == "Mo" && $2 == co && $3 == T) print $4}' numbers.txt); # these three are for the EAM potential
        # nMo=$(awk -v co=${i} -v T=${j} '{if ($1 == "Mo" && $2 == co && $3 == T) print $5}' numbers.txt);
        nXe=0
        n=$((${nU}+${nMo}+${nXe}))
        awk -v nU=${nU} 'NR == 1 {print $0} NR == 2 {print $0,$NF"*"nU} NR > 2 {print $0" "$NF*nU}' MSD_U_T${j}_${i}Mo.txt > tmp_${i}_${j}_U.txt;
        awk -v nMo=${nMo} 'NR == 1 {print $0} NR == 2 {print $0,$NF"*"nMo} NR >2 {print $0" "$NF*nMo}' MSD_${i}Mo_T${j}.txt > tmp_${i}_${j}_Mo.txt;
        awk -v n=${n} 'NR == 1 {print $0} NR == 2 {print $0,$NF"*"n} NR > 2 {print $0" "$NF*n}' MSD_all_T${j}_${i}Mo.txt > tmp_${i}_${j}_all.txt;
        mv tmp_${i}_${j}_U.txt MSD_U_T${j}_${i}Mo.txt;
        mv tmp_${i}_${j}_Mo.txt MSD_${i}Mo_T${j}.txt;
        mv tmp_${i}_${j}_all.txt MSD_all_T${j}_${i}Mo.txt
      # elif [[ $(qalc -t ${i} == 0) -eq 1 ]] && [[ $(qalc -t ${k} != 0) -eq 1 ]]; then
      #   nU=$(awk -v co=${k} -v T=${j} '{if ($1 == "Xe" && $2 == co && $3 == T) print $4}' numbers.txt);
      #   nXe=$(awk -v co=${k} -v T=${j} '{if ($1 == "Xe" && $2 == co && $3 == T) print $5}' numbers.txt);
      #   nMo=0
      #   n=$((${nU}+${nMo}+${nXe}))
      #   awk -v nU=${nU} 'NR == 1 {print $0} NR == 2 {print $0,$NF"*"nU} NR > 2 {print $0" "$NF*nU}' MSD_U_T${j}_${k}Xe.txt > tmp_${k}_${j}_U.txt;
      #   awk -v nXe=${nXe} 'NR == 1 {print $0} NR == 2 {print $0,$NF"*"nXe} NR >2 {print $0" "$NF*nXe}' MSD_${k}Xe_T${j}.txt > tmp_${k}_${j}_Xe.txt;
      #   awk -v n=${n} 'NR == 1 {print $0} NR == 2 {print $0,$NF"*"n} NR > 2 {print $0" "$NF*n}' MSD_all_T${j}_${k}Xe.txt > tmp_${k}_${j}_all.txt;
      #   mv tmp_${k}_${j}_U.txt MSD_U_T${j}_${k}Xe.txt;
      #   mv tmp_${k}_${j}_Xe.txt MSD_${k}Xe_T${j}.txt;
      #   mv tmp_${k}_${j}_all.txt MSD_all_T${j}_${k}Xe.txt
      # else
      #   continue
      # fi;
    done
  # done;
done

r="/media/jarinf/Research Backup/Research1_backup/U/moly_effect/grain_growth/gamma/110/ternary_eam"
cd "${r}"
# for i in {0.{25,50,75},1}; do
for i in 3; do # {1,3,5,10,14,17,21,30,50}
  for j in {20,30,45}; do
    for k in {1300..1400..50}; do
      for l in large_r/dir_test_{1..3}; do
        if [[ -d "${i}at%/${j}degree/T${k}/${l}" ]]; then
          (
          cd "${i}at%/${j}degree/T${k}/${l}"
          echo "$(ls *.dump | wc -l) ${j} 2 1.207 0.0 3.463 bcc" > find_grains_input.txt
          echo -e "0 0 1\n1 -1 0\n1 1 0" >> find_grains_input.txt
          ls -v *.dump >> find_grains_input.txt

          echo "data.txt ${k} $(qalc -t ${i}/100) $(awk 'NR==7 {print $2}' "${r}/${i}at%/${j}degree/LAMMPS_U${i}Mo_110_pTernaryEAM_${j}degree_r100A_annealed.dat") 3.463 bcc" > grain_area_input.txt
          # echo "data.txt ${k} $(qalc -t ${i}/100) $(awk 'NR==7 {print $2}' "${r}/${i}at%/${j}degree/LAMMPS_U${i}Xe_110_pTernaryEAM_${j}degree_r100A_annealed.dat") 3.463 bcc" > grain_area_input.txt
          find_grains find_grains_input.txt -e 10
          calculate_grain_area grain_area_input.txt -p 6
          )
        fi
      done
    done
  done
done

reg="_U(([0-9]*)Mo)?(([0-9]*(\.([0-9]*)))?Xe)?_"
for i in "${anneals_list[@]}"; do
  SEED=$(bash -c 'echo "${RANDOM}"')
  if [[ ${i} =~ ${reg} ]]; then
    ii=1
    n=${#BASH_REMATCH[*]}
    while [[ "${ii}" -lt $n ]]; do
      echo "  capture[$ii]: ${BASH_REMATCH[$ii]}"
      let ii++
    done
  fi
  vars_list="-var SEED ${SEED} -var base ${i/%.dat/} -in anneal.in"
  mpirun ${LAMMPS_EXEC} ${vars_list} > output_${i/%.dat/}.txt
done

reg="_U(([0-9]*)Mo)?(([0-9]*(\.([0-9]*)))?Xe)?_"
for i in *.dat; do
  if [[ "${i}" =~ ${reg} ]]; then
    ii=1
    n="${#BASH_REMATCH[*]}"
    while [[ "${ii}" -lt "$n" ]]; do
      if [[ -z "${BASH_REMATCH[${ii}]}" ]]; then
        let ii++
        continue
      fi
      echo "  capture[${ii}]: ${BASH_REMATCH[${ii}]}"
      let ii++
    done
  fi
done

for i in $(fd dir_test); do
  (
  cd ${i}
  cx=$(qalc -t $(awk 'NR == 5 {print $2}' ../../../LAMMPS_*)/2)
  cy=$(qalc -t $(awk 'NR == 6 {print $2}' ../../../LAMMPS_*)/2)
  concentration_analysis $(ls -v1 *.dump) -t 2 3 -s 3.463 -d r -c ${cx} ${cy}
  )
done

for i in $(fd find_grains_input.txt -E dir_[0-9]); do
  a=$(awk 'NR==1 {print $6}' ${i})
  if [[ $(qalc -t "${a} != 3.542") ]]; then
    (
    cd $(dirname ${i})
    sed -i '1s/3.463/3.542/' find_grains_input.txt
    sed -i '1s/3.463/3.542/' grain_area_input.txt
    mv area_data.txt area_data_wrong.txt
    mv data.txt data_wrong.txt
    find_grains find_grains_input.txt -e 0
    calculate_grain_area grain_area_input.txt -p 3
    mobility_plots.py area_data.txt -F
    )
  fi
done


#Useful for compressing the the impure U grain growth simulation data tp the NAS
for i in U[0-9]*.tar.bz2; do
  axis=$(echo "${i}" | awk -F'_' '{print $2}')
  mis=$(echo "${i}" | awk -F'_' '{print $3}')
  T=$(echo "${i}" | awk -F'_' '{print $4}')
  syst=$(echo "${i}" | awk -F'_' '{print $1}')
  rgx="U([0-9]+(.[0-9][0-9])?(Mo|Xe))([0-9].[0-9][0-9]Xe)?"

  if [[ ${syst} =~ ${rgx} ]]; then
    if [[ ${BASH_REMATCH[3]} == "Mo" ]] && [[ -z ${BASH_REMATCH[4]} ]]; then
      sys_dir="moly_effect"
      cMo=$(echo "${BASH_REMATCH[1]}" | cut -d "M" -f 1)
      cXe=0
      imp_dir="${cMo}at%"
    elif [[ ${BASH_REMATCH[3]} == "Mo" ]] && [[ ! -z ${BASH_REMATCH[4]} ]]; then
      sys_dir="moly_xenon_effect"
      cMo=$(echo "${BASH_REMATCH[1]}" | cut -d "M" -f 1)
      cXe=$(echo "${BASH_REMATCH[4]}" | cut -d "X" -f 1)
      imp_dir="${cMo}at%Mo/${cXe}at%Xe"
    else
      sys_dir="xenon_effect"
      cMo=0
      cXe=$(echo "${BASH_REMATCH[1]}" | cut -d "X" -f 1)
      imp_dir="${cXe}at%"
    fi
  else
    echo "Regex matching failed for ${i}"
    continue
  fi

  archive="${BASH_REMATCH[0]}_${axis}_${mis}_ternary_eam_${T}_large_r_results.tar.7z"
  if [[ -f /nfs/home/Research/U/impurity_effect/grain_growth/${archive} ]]; then
    continue
  fi
  tar -xjvf "${i}" --strip-components=3 # extract the tar.bz2 file
  (
  cd "U/${sys_dir}/grain_growth/gamma/${axis}/ternary_eam/${imp_dir}/${mis}/${T}/large_r" || continue
  for j in dir_*; do
    (
    cd "${j}" || continue
    if [[ -f snapshots.tar.bz2 ]] &&  [[ ! -f snapshots.7z ]]; then
      tar -xjvf snapshots.tar.bz2 # extract the tar.bz2 file containing the dump files
      7za snapshots.7z ./*.dump && rm snapshots.tar.bz2 ./*.dump # we don't want duplicated data in the resulting archive, so remove the old archive as well as the (now compressed) files
    elif [[ -f snapshots.7z ]] && [[ -f snapshots.tar.bz2 ]]; then
      rm snapshots.tar.bz2 *.dump
    fi
    )
  done
  tar -cvf - $(fd . -t f -E '*_interface.dat' -E mobility_plot.png -E '*.dump') | 7za -si "/nfs/home/Research/U/impurity_effect/grain_growth/${archive}"
  )
done

# For plotting the velocity vs the force
for i in {100,110,111,112}/{{20,45}degree,sigma7}/{Basak,Cooper}/T{{20..33}00,3050}/{large,medium}_r/dir_final_{1..5}; do
  if [[ ! -d ${i} ]]; then
    continue
  fi
  T=$(echo ${i} | awk -F'/' '{print $4}' | cut -c 2-)
  rdir=$(echo ${i} | awk -F'/' '{print $5}')
  if [[ ${rdir} == "large_r" ]]; then
    r=100
  elif [[ ${rdir} == "medium_r" ]]; then
    r=75
  else
    echo "Unknown rdir: ${rdir}"
    continue
  fi
  (
  cd ${i}
  if [[ ! -f "force_velocity_data.txt" ]]; then
    continue
  fi
  gnuplot -e T=${T} -e "el='UO2'" -e r=${r} plot_vel_vs_force.plt
  )
done

# For this specific case, we draw four (resized) images on top of an original image, and also draw a dashed circle
convert -gravity center Basak_111_20degree_T2800_dir1_vector_plot.png  \( cw_arrows.png -resize 50x50 \) -geometry -330+10 -composite  \( cw_arrows.png -resize 50x50 \) -geometry +320+10 -composite  \( cw_arrows.png -resize 50x50 \) -geometry +18+315 -composite  \( cw_arrows.png -resize 50x50 \) -geometry +10-325 -composite -fill none -stroke black -gravity center -draw 'translate 450,370 stroke-dasharray 5 5 ellipse 100,100 300,300 0,360'  Basak_111_20degree_T2800_dir1_vector_plot_with_indicator.png

for i in T{17..33}00/large_r/dir_{1..6}/; do
  if [[ -d "${i}" ]] && [[ ! -f "${i}/MSD_U.dat" ]]; then
    (
    cd ${i}
    extract snapshots.7z
    python ~/projects/scripts/ovito/calculate_MSD.py && rm *.dump
    )
  fi
done

for msd in xy z; do
  for j in {1..5}; do
    for i in {20..30}00 3050; do
      if [[ ! -d T${i}/large_r/dir_final_${j} ]];
      then continue;
      fi;
      (
      cd T${i}/large_r/dir_final_${j};
      python ~/projects/scripts/python/calculate_diffusion.py --show ${msd} | awk '{print $5,$9}');
    done;
    echo " ";
  done;
  echo "-------------------------------------------------------";
done

for i in $(fd dir_ -t d); do
  (
  cd ${i}
  if [[ ! -f "rdf_initial.txt" ]]; then
    python ~/projects/scripts/ovito/calculate_RDF.py
  fi
  if [[ ! -f "cluster_results_type2_initial.txt" ]]; then
    python ~/projects/scripts/ovito/cluster_analysis.py bcc 3.542
  fi
  if [[ ! -f "force_velocity_data.txt" ]]; then
    r1="T(1[0-4][05]0)"
    r2="(([0-9](.[025][05])?)at%(Mo|Xe)?)"
    if [[ "${i}" =~ ${r1} ]]; then
      T=${BASH_REMATCH[1]}
    else
      echo "T not found for ${i}"
    fi
    if [[ "${i}" =~ ${r2} ]]; then
      c=$(echo "scale=2;${BASH_REMATCH[2]}/100" | bc)
    else
      echo "c not found for ${i}"
    fi
    if [[ -f "log.lammps" ]]; then
      l=$(grep "orthogonal box" log.lammps | awk '{print $6,$10}' | awk -F')' '{print $2-$1}')
    elif compgen -G "output*.txt" > /dev/null; then
      l=$(grep "orthogonal box" log.lammps | awk '{print $6,$10}' | awk -F')' '{print $2-$1}')
    else
      echo "Unable to determine box size (missing log.lammps and output*.txt)"
    fi
    p=$(pwd | awk -F'/' '{print $7}' | sed 's/_effect//')
    if [[ ${p} == "xenon" ]]; then
      pnum=6
      calculate_force_and_velocity.py ${T} ${c} ${l} 3.542 1.0 -p ${pnum} -u 2
    elif [[ ${p} == "moly" ]] || [[ ${p} == "moly_xenon" ]]; then
      pnum=2
      calculate_force_and_velocity.py ${T} ${c} ${l} 3.542 1.0 -p ${pnum} -u 2
    else
      echo "Unrecognized potential ${p}"
    fi
  fi
  if ! compgen -G "*_displacement_data.txt" > /dev/null; then
    f=0.dump
    s=$(ls -v *.dump | head -n "$(qalc -t "round(0.75*$(tail -n1 slope_calc.txt | awk -F':' '{print $2}'))")" | tail -n1)
    track_atoms 0.dump --xlo 0 --xhi 1 --ylo 0 --yhi 1 --zlo 0 --zhi 1 -i 1 --print-atom-ids > tmp
    calculate_displacement ${f} ${s} --ids-file tmp
    rm tmp
  fi
  if [[ ! -f "trajectories_type2.png" ]]; then
    python3.10 ~/projects/scripts/ovito/generate_trajectory_lines.py
  fi
  if [[ ! -f "diffusion_results.txt" ]]; then
    python ~/projects/scripts/python/calculate_diffusion.py $(ls MSD_*.dat | tail -n +2) --post-growth-diff
  fi
  )
done
