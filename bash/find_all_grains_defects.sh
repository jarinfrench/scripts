#! /bin/bash

trap ctrl_c INT

function ctrl_c() {
  # Force exiting of the script, rather than just the current executable
  exit 1
}

lattice_param=5.453
cut100=0.866
cut110=0.866
cut111=1.11237
cutoff=1.2
gamma=1.6


cd /media/jarinf/Seagate\ Expansion\ Drive/uo2/defect_evolution/110/45degree/
#
# for i in 110 111; do
#   cd ${i}/45degree
  for j in $(ls -vd T*/); do
    cd ${j}
    for i in large_r small_r tiny_r; do
      if [ -d ${i} ]; then
        cd ${i}
      else
        continue
      fi
    for k in $(ls -vd dir_*/); do
      cd ${k};
      if [ -d "interfaces" ]; then
        echo "Directory \"interfaces\" already found in ${j}${i}/${k}"
        if ! [ -d "tracked_data" ]; then
          mkdir tracked_data
          track_atoms `ls -v *.dump`
          mv *_tracked* tracked_data
        fi
        t=$(echo $j | cut -c 2- | awk -F '/' '{print $1}')
        height=$(cat UO2_minimized_* | head -n 7 | tail -n 1 | awk '{print $2}')
        potential=2
        echo "data.txt ${t} ${height} ${lattice_param} fcc" > grain_area_input.txt
        calculate_grain_area grain_area_input.txt -p ${potential}
        # if ! [ -a "force_velocity_data.txt" ]; then
        #   calculate_force_and_velocity.py $t ${height} ${lattice_param} ${gamma} -p $potential -g
        # fi
        # find_new_positions 0.dump $(ls -v *.dump | tail -n 1) tracked_positions.dat
        echo -e "\n"
      else
        if [ "$(ls *.dump 2>/dev/null | wc -l)" -eq 0 ]; then
          cd ..
          continue
        fi
        # cutoff=1.25

        # if [ "$i" -eq 100 ]; then
        #   echo "$(ls *.dump | wc -l) 45.00 3 ${cut100} ${cutoff} ${lattice_param}" > base_vals.txt
        #   echo "1 0 0" >> base_vals.txt
        #   echo "0 1 0" >> base_vals.txt
        #   echo "0 0 1" >> base_vals.txt
        # elif [ "$i" -eq 110 ]; then
          echo "$(ls *.dump | wc -l) 45.00 3 ${cut110} ${cutoff} ${lattice_param}" > base_vals.txt
          echo "0 0 1" >> base_vals.txt
          echo "1 -1 0" >> base_vals.txt
          echo "1 1 0" >> base_vals.txt
        # elif [ "$i" -eq 111 ]; then
        #   echo "$(ls *.dump | wc -l) 38.20 3 ${cut111} ${cutoff} ${lattice_param}" > base_vals.txt
        #   echo "1 -1 0" >> base_vals.txt
        #   echo "1 1 -2" >> base_vals.txt
        #   echo "1 1 1" >> base_vals.txt
        # fi

        ls -v *.dump >> base_vals.txt
        find_grains base_vals.txt -i 2
        mkdir interfaces
        mv *_interface* interfaces

        t=$(echo ${j} | cut -c 2- | awk -F '/' '{print $1}')
        height=$(cat UO2_minimized_* | head -n 7 | tail -n 1 | awk '{print $2}')
        potential=2
        echo "data.txt ${t} ${height} ${lattice_param} fcc" > grain_area_input.txt
        calculate_grain_area grain_area_input.txt -p ${potential}
        echo -e "\n"
      fi
      cd ..
    done
    cd ../
    done
    cd ../
  done


#   cd ../..
# done
