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
            for m in $(ls -d dir*/); do
              cd ${m}
              if [ -d "interfaces" ]; then
                echo "Directory \"interfaces\" already found in ${i}/45degree/${j}/${k}${l}/${m}"
                if ! [ -a "area_data.txt" ]; then
                  t=$(echo $k | cut -c 2-)
                  height=$(cat UO2_minimized_* | head -n 7 | tail -n 1 | awk '{print $2}')
                  if [ "$j" == "basak" ]; then
                    potential=2
                  elif [ "$j" == "eam" ]; then
                    potential=1
                  fi
                  calculate_grain_area data.txt $t $height 5.453 $potential
                else
                  echo -e "\033[0;32m\tArea data file already found in ${i}/45degree/${j}/${k}${l}/${m}\033[0m"
                fi
              else
                if [ "$i" -eq 100 ]; then
                  echo "$(ls *.dump | wc -l) 45.00 2 0.866 5.453" > base_vals.txt
                  echo "1 0 0" >> base_vals.txt
                  echo "0 1 0" >> base_vals.txt
                  echo "0 0 1" >> base_vals.txt
                elif [ "$i" -eq 110 ]; then
                  echo "$(ls *.dump | wc -l) 45.00 2 0.866 5.453" > base_vals.txt
                  echo "0 0 1" >> base_vals.txt
                  echo "1 -1 0" >> base_vals.txt
                  echo "1 1 0" >> base_vals.txt
                elif [ "$i" -eq 111 ]; then
                  echo "$(ls *.dump | wc -l) 45.00 2 1.2 5.453" > base_vals.txt
                  echo "1 -1 0" >> base_vals.txt
                  echo "1 1 -2" >> base_vals.txt
                  echo "1 1 1" >> base_vals.txt
                fi
                ls -v *.dump >> base_vals.txt
                find_grains base_vals.txt
                mkdir interfaces
                mv *_interface* interfaces
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
