#! /bin/bash

trap ctrl_c INT

function ctrl_c() {
  # Force exiting of the script, rather than just the current executable
  exit 1
}

cd /media/jarinf/Seagate\ Expansion\ Drive/U/xenon_effect/gamma/


for m in 100; do # 110 111; do
  cd ${m}
  for i in ternary_eam; do # adp; do
    cd ${i}
    for j in 30degree 45degree; do # sigma7; do
      if [[ -d ${j} ]]; then
        cd ${j}
      else
        continue
      fi
      for z in 0.5 1.0 1.5 2.0 2.5 5.0 7.5 10.0; do
        cd ${z}%
        for k in $(ls -vd T*/); do
          cd ${k}large_r
          for l in $(ls -vd dir*/); do
            cd ${l}
            if [[ "${i}" =~ ^adp$ ]]; then
              lattice_param=3.52
              potential=6
              rcut=0.933
              fi_cut=1.4
            elif [[ "${i}" =~ ^ternary_eam$ ]]; then
              lattice_param=3.542
              potential=9
              rcut=0.933
              fi_cut=1.4
            fi

            echo -e "\nProcessing directory $(pwd)"

            if [ -d "interfaces" ]; then
              echo "Directory \"interfaces\" already found"
              if [[ $(ls *.dump | wc -l) -gt 0 ]]; then
                7z a U_${z}Xe_${m}_${i}_${k%/}_large_r_${l%/}_dump.7z *.dump
                if [ $? -eq 0 ]; then
                  rm *.dump
                else
                  RED='\033[0;31m'
                  echo -e "${RED}Error creating zip file"
                fi
              fi
              t=$(echo ${k} | cut -c 2- | awk -F '/' '{print $1}')
              height=$(cat U_minimized_* | head -n 7 | tail -n 1 | awk '{print $2}')
              echo "data.txt ${t} ${height} ${lattice_param} bcc" > grain_area_input.txt
              calculate_grain_area grain_area_input.txt -p ${potential}
              cd ..
            else
              if [[ "${m}" -eq 100 ]]; then
                if [[ "${j}" =~ ^30degree$ ]]; then
                  theta=7.00
                elif [[ "${j}" =~ ^45degree$ ]]; then
                  theta=0.00
              #   elif [[ "${j}" =~ ^sigma7$ ]]; then
              #     theta=0.00
                else
                  echo "Error: Unrecognized directory"
                  continue
                fi
                echo "$(ls *.dump | wc -l) ${theta} 2 ${rcut} ${fi_cut} ${lattice_param}" > base_vals.txt
                echo "1 0 0" >> base_vals.txt
                echo "0 1 0" >> base_vals.txt
                echo "0 0 1" >> base_vals.txt
              # elif [[ "${m}" -eq 110 ]]; then
              #   if [[ "${j}" =~ ^30degree$ ]]; then
              #     theta=13.00
              #   elif [[ "${j}" =~ ^45degree$ ]]; then
              #     theta=8.00
              #   elif [[ "${j}" =~ ^sigma7$ ]]; then
              #     theta=0.00
              #   else
              #     echo "Error: Unrecognized directory"
              #     continue
              #   fi
              #   echo "$(ls *.dump | wc -l) ${theta} 1 ${rcut} ${fi_cut} ${lattice_param}" > base_vals.txt
              #   echo "0 0 1" >> base_vals.txt
              #   echo "1 -1 0" >> base_vals.txt
              #   echo "1 1 0" >> base_vals.txt
              # elif [[ "${m}" -eq 111 ]]; then
              #   if [[ "${j}" =~ ^30degree$ ]]; then
              #     theta=77.00
              #   elif [[ "${j}" =~ ^45degree$ ]]; then
              #     theta=74.20
              #   elif [[ "${j}" =~ ^sigma7$ ]]; then
              #     theta=75.60
              #   else
              #     echo "Error: Unrecognized directory"
              #     continue
              #   fi
              #   echo "$(ls *.dump | wc -l) ${theta} 1 ${rcut} ${fi_cut} ${lattice_param}" > base_vals.txt
              #   echo "1 -1 0" >> base_vals.txt
              #   echo "1 1 -2" >> base_vals.txt
              #   echo "1 1 1" >> base_vals.txt
              fi
              ls -v *.dump >> base_vals.txt

              find_grains base_vals.txt -e 10 -q # interface files written every 10 files, and suppress warnings.
              mkdir interfaces
              mv *_interface.dat interfaces

              if [[ $(ls *.dump | wc -l) -gt 0 ]]; then
                7z a U_${z}Xe_${m}_${i}_${j}_${k%/}_large_r_${l%/}_dump.7z *.dump
                if [ $? -eq 0 ]; then
                  rm *.dump
                else
                  RED='\033[0;31m'
                  echo -e "${RED}Error creating zip file"
                fi
              fi


              # if ! [ -d "tracked_data" ]; then
              #   mkdir tracked_data
              #   track_atoms 10000.dump $(ls -v *.dump | tail -n +3)
              #   mv *_tracked.dat tracked_data
              # fi

              t=$(echo ${k} | cut -c 2- | awk -F '/' '{print $1}')
              height=$(cat U_minimized_* | head -n 7 | tail -n 1 | awk '{print $2}')
              echo "data.txt ${t} ${height} ${lattice_param} bcc" > grain_area_input.txt
              calculate_grain_area grain_area_input.txt -p ${potential}
              cd ..
            fi
          done
          cd ../../
        done
        cd ../
      done
      cd ../
    done
    cd ../
  done
  cd ../
done

        # if [ -d "interfaces" ]; then
        #   echo "Directory \"interfaces\" already found in ${i}/${j}${k}"
        #
        #
        #   if ! [ -a "force_velocity_data.txt" ]; then
        #     calculate_force_and_velocity.py $t ${height} ${lattice_param} ${gamma} -p $potential -g
        #   fi
