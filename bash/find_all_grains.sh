#! /bin/bash

cd /media/jarinf/Seagate\ Expansion\ Drive/UO2/grain_growth

for i in 100 110 111; do
  cd ${i}/45degree
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
                calculate_grain_area data.txt $t ${height} 5.453 ${potential}
                calculate_force_and_velocity.py $t ${height} 5.453 -p ${potential}
                echo -e "\n"
              else
                if [ "$(ls *.dump 2>/dev/null | wc -l)" -eq 0 ]; then
                  cd ..
                  continue
                fi
                if [ "$i" -eq 100 ]; then
                  echo "$(ls *.dump | wc -l) 45.00 2 0.866 1.25 5.453" > base_vals.txt
                  echo "1 0 0" >> base_vals.txt
                  echo "0 1 0" >> base_vals.txt
                  echo "0 0 1" >> base_vals.txt
                elif [ "$i" -eq 110 ]; then
                  echo "$(ls *.dump | wc -l) 45.00 2 0.866 1.25 5.453" > base_vals.txt
                  echo "0 0 1" >> base_vals.txt
                  echo "1 -1 0" >> base_vals.txt
                  echo "1 1 0" >> base_vals.txt
                elif [ "$i" -eq 111 ]; then
                  echo "$(ls *.dump | wc -l) 45.00 2 1.2 1.25 5.453" > base_vals.txt
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
                calculate_grain_area data.txt $t ${height} 5.453 ${potential}
                calculate_force_and_velocity.py $t ${height} 5.453 -p ${potential}
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
  cd ../..
done
