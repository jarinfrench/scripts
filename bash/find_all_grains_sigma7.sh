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


cd /media/jarinf/Seagate\ Expansion\ Drive/uo2/grain_growth/111/sigma7




for k in $(ls -vd T*/); do
  cd ${k}large_r
  for m in $(ls -d dir*/); do
    cd ${m}
    if [ -d "interfaces" ]; then
      echo "Directory \"interfaces\" already found in ${k}large_r/${m}"
      if ! [ -d "tracked_data" ]; then
        ls -v *.dump >> tracked_input.txt
        mkdir tracked_data
        find_new_positions tracked_input.txt
        mv *_tracked* tracked_data
      fi
      t=$(echo $k | cut -c 2- | awk -F '/' '{print $1}')
      height=$(cat UO2_minimized_* | head -n 7 | tail -n 1 | awk '{print $2}')
      potential=2
      echo "data.txt ${t} ${height} ${lattice_param} fcc" > grain_area_input.txt
      calculate_grain_area grain_area_input.txt -p ${potential}
      # if ! [ -a "force_velocity_data.txt" ]; then
      #   calculate_force_and_velocity.py $t ${height} ${lattice_param} ${gamma} -p $potential -g
      # fi
      find_new_positions $(ls -v *.dump | head -n 3 | tail -n 1) $(ls -v *.dump | tail -n 1) tracked_positions.dat
      echo -e "\n"
    else
      if [ "$(ls *.dump 2>/dev/null | wc -l)" -eq 0 ]; then
        cd ..
        continue
      fi
      cutoff=1.25

      echo "$(ls *.dump | wc -l) 20.00 2 ${cut111} ${cutoff} ${lattice_param}" > base_vals.txt
      echo "1 -1 0" >> base_vals.txt
      echo "1 1 -2" >> base_vals.txt
      echo "1 1 1" >> base_vals.txt

      ls -v *.dump >> base_vals.txt
      ls -v *.dump >> tracked_input.txt
      find_grains base_vals.txt -i 2
      find_new_positions tracked_input.txt
      mkdir interfaces tracked_data
      mv *_interface* interfaces
      mv *_tracked* tracked_data

      t=$(echo $k | cut -c 2- | awk -F '/' '{print $1}')
      height=$(cat UO2_minimized_* | head -n 7 | tail -n 1 | awk '{print $2}')
      potential=2
      echo "data.txt ${t} ${height} ${lattice_param} fcc" > grain_area_input.txt
      calculate_grain_area grain_area_input.txt -p ${potential}
      find_new_positions 0.dump $(ls -v *.dump | tail -n 1) tracked_positions.dat
      echo -e "\n"
    fi
    cd ..
  done
  cd ../..
done
