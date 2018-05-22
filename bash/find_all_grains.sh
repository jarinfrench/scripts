#! /bin/bash

trap ctrl_c INT

function ctrl_c() {
  # Force exiting of the script, rather than just the current executable
  exit 1
}

read -p "Please enter the element: " element
if [[ "${element,,}" =~ ^(al|cu|fe)$ ]]; then
  element=${element,,}
  element=${element^}
elif [[ "${element,,}" =~ ^uo2$ ]]; then
  element=${element^^}
else
  echo "$element is not recognized.  Please enter either UO2, Fe, Cu, or Al"
  exit 3
fi

if [[ "${element}" =~ ^Fe$ ]]; then
  potential=5
  lattice_param=2.85532
  cut100=0.933
  cut110=0.933
  cut111=0.933
  cutoff=1.6
  gamma=1.0
elif [[ "${element}" =~ ^Al$ ]]; then
  potential=4
  lattice_param=4.032
  cut100=0.866
  cut110=0.866
  cut111=1.2
  cutoff=1.25
  gamma=0.5
elif [[ "${element}" =~ ^Cu$ ]]; then
  potential=3
  lattice_param=3.615
  cut100=0.866
  cut110=0.866
  cut111=1.2
  cutoff=1.25
  gamma=0.9
elif [[ "${element}" =~ ^UO2$ ]]; then
  lattice_param=5.453
  cut100=0.866
  cut110=0.866
  cut111=1.2
  cutoff=1.25
  gamma=1.6
else
  potential=0
  read -p "Please enter the lattice parameter of the material: " lattice_param
  read -p "Please enter the random grain boundary energy of this material: " gamma
  while true; do
    read -p "Is the material FCC (1), or BCC (2)? " crystal_type
    if [[ "${crystal_type}" == 1 ]] || [[ "${crystal_type}" == 2 ]]; then
      break
    fi
  done
  if [[ "${crystal_type}" == 1 ]]; then
    cut100=0.866
    cut110=0.866
    cut111=1.2
    cutoff=1.25
  elif [[ "${crystal_type}" == 2 ]]; then
    cut100=0.933
    cut110=0.933
    cut111=0.933
    cutoff=1.6
  else
    echo "Error determining crystal type."
    exit 4
  fi
fi

cd /media/jarinf/Seagate\ Expansion\ Drive/${element}/grain_growth

for i in 100 110 111; do
  cd ${i}
  if ! [[ "${element}" =~ ^UO2$ ]]; then
    for j in $(ls -vd T*/); do
      cd ${j}
      for k in $(ls -d dir*/); do
        cd ${k}
        echo -e "\nProcessing directory ${i}/${j}${k}"
        if [ -d "interfaces" ]; then
          echo "Directory \"interfaces\" already found in ${i}/${j}${k}"
          t=$(echo $j | cut -c 2- | awk -F '/' '{print $1}')
          height=$(cat ${element}_minimized_* | head -n 7 | tail -n 1 | awk '{print $2}')
          calculate_grain_area data.txt $t $height ${lattice_param} $potential
          if ! [ -a "force_velocity_data.txt" ]; then
            calculate_force_and_velocity.py $t ${height} ${lattice_param} ${gamma} -p $potential -g
          fi
          find_new_positions 0.dump $(ls -v *.dump | tail -n 1) tracked_positions.dat
        else
          if [ "$(ls *.dump 2>/dev/null | wc -l)" -eq 0 ]; then
            cd ..
            continue
          fi
          if [ "$i" -eq 100 ]; then
            echo "$(ls *.dump | wc -l) 45.00 1 ${cut100} ${cutoff} ${lattice_param}" > base_vals.txt
            echo "1 0 0" >> base_vals.txt
            echo "0 1 0" >> base_vals.txt
            echo "0 0 1" >> base_vals.txt
          elif [ "$i" -eq 110 ]; then
            echo "$(ls *.dump | wc -l) 45.00 1 ${cut110} ${cutoff} ${lattice_param}" > base_vals.txt
            echo "0 0 1" >> base_vals.txt
            echo "1 -1 0" >> base_vals.txt
            echo "1 1 0" >> base_vals.txt
          elif [ "$i" -eq 111 ]; then
            echo "$(ls *.dump | wc -l) 45.00 1 ${cut111} ${cutoff} ${lattice_param}" > base_vals.txt
            echo "1 -1 0" >> base_vals.txt
            echo "1 1 -2" >> base_vals.txt
            echo "1 1 1" >> base_vals.txt
          fi
          ls -v *.dump >> base_vals.txt
          find_grains base_vals.txt
          mkdir interfaces
          mv *_interface* interfaces
          if ! [ -a "area_data.txt" ]; then
            t=$(echo $j | cut -c 2- | awk -F '/' '{print $1}')
            height=$(cat ${element}_minimized_* | head -n 7 | tail -n 1 | awk '{print $2}')
            calculate_grain_area data.txt $t ${height} ${lattice_param} $potential
            if ! [ -a "force_velocity_data.txt" ]; then
              calculate_force_and_velocity.py $t ${height} ${lattice_param} ${gamma} -p $potential -g
            fi
            find_new_positions 0.dump $(ls -v *.dump | tail -n 1) tracked_positions.dat
          else
            echo -e "\033[0;32m\tArea data file already found in ${i}/${j}${k}\033[0m"
          fi
        fi
        cd ..
      done
      cd ..
    done
    cd ..
  else
    cd 45degree
    for j in basak eam; do
      cd ${j}
      for k in $(ls -vd T*/); do
        cd ${k}
        for l in large_r medium_r small_r tiny_r; do
          if [ -d "$l" ]; then
            cd ${l}
            #if ! [ -d "dir*" ]; then
            #  continue
            #fi
            for m in $(ls -d dir*/ 2>/dev/null); do
              cd ${m}
              if [ -d "interfaces" ]; then
                echo "Directory \"interfaces\" already found in ${i}/45degree/${j}/${k}${l}/${m}"
                t=$(echo $k | cut -c 2- | awk -F '/' '{print $1}')
                height=$(cat UO2_minimized_* | head -n 7 | tail -n 1 | awk '{print $2}')
                if [ "$j" == "basak" ]; then
                  potential=2
                elif [ "$j" == "eam" ]; then
                  potential=1
                fi
                calculate_grain_area data.txt $t ${height} ${lattice_param} ${potential}
                if ! [ -a "force_velocity_data.txt" ]; then
                  calculate_force_and_velocity.py $t ${height} ${lattice_param} ${gamma} -p $potential -g
                fi
                find_new_positions 0.dump $(ls -v *.dump | tail -n 1) tracked_positions.dat
                echo -e "\n"
              else
                if [ "$(ls *.dump 2>/dev/null | wc -l)" -eq 0 ]; then
                  cd ..
                  continue
                fi
                if [ "$i" -eq 100 ]; then
                  echo "$(ls *.dump | wc -l) 45.00 2 ${cut100} ${cutoff} ${lattice_param}" > base_vals.txt
                  echo "1 0 0" >> base_vals.txt
                  echo "0 1 0" >> base_vals.txt
                  echo "0 0 1" >> base_vals.txt
                elif [ "$i" -eq 110 ]; then
                  echo "$(ls *.dump | wc -l) 45.00 2 ${cut110} ${cutoff} ${lattice_param}" > base_vals.txt
                  echo "0 0 1" >> base_vals.txt
                  echo "1 -1 0" >> base_vals.txt
                  echo "1 1 0" >> base_vals.txt
                elif [ "$i" -eq 111 ]; then
                  echo "$(ls *.dump | wc -l) 45.00 2 ${cut111} ${cutoff} ${lattice_param}" > base_vals.txt
                  echo "1 -1 0" >> base_vals.txt
                  echo "1 1 -2" >> base_vals.txt
                  echo "1 1 1" >> base_vals.txt
                fi
                ls -v *.dump >> base_vals.txt
                find_grains base_vals.txt
                mkdir interfaces
                mv *_interface* interfaces

                t=$(echo $k | cut -c 2- | awk -F '/' '{print $1}')
                height=$(cat UO2_minimized_* | head -n 7 | tail -n 1 | awk '{print $2}')
                if [ "$j" == "basak" ]; then
                  potential=2
                elif [ "$j" == "eam" ]; then
                  potential=1
                fi
                calculate_grain_area data.txt $t ${height} ${lattice_param} ${potential}
                if ! [ -a "force_velocity_data.txt" ]; then
                  calculate_force_and_velocity.py $t ${height} ${lattice_param} ${gamma} -p $potential -g
                fi
                find_new_positions 0.dump $(ls -v *.dump | tail -n 1) tracked_positions.dat
                echo -e "\n"
              fi
              cd ..
            done
            cd ..
          fi
        done
        cd ..
      done
      cd ..
    done
    cd ..
  cd ..
  fi
done
