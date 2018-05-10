#! /bin/bash

cd /media/jarinf/Seagate\ Expansion\ Drive/UO2/grain_growth

for i in 100 110 111; do
  cd ${i}/45degree
  for j in basak eam; do
    cd ${j}
    # echo "${i}/45degree/${j}: $temps"
    for k in $(ls -vd T*/); do
      cd ${k}
        for l in large_r medium_r small_r tiny_r; do
          if [ -d "$l" ]; then
            cd ${l}
            for m in $(ls -d dir*/); do
              cd ${m}
              if [ -d "interfaces" ]; then
                cd interfaces
              else
                echo -e "\033[0;31mUnable to find directory \"interfaces\" in directory ${i}/45degree/${j}/${k}${l}/$m\033[0m"
                cd ..
                continue
              fi
              if [ -d "displacement_data" ]; then
                echo "directory \"displacement_data\" already found in directory ${i}/45degree/${j}/${k}${l}/${m}interfaces"
                cd ../..
                continue
              fi
              ls -v *_interface.dat > displacement_input.txt
              calculate_displacement displacement_input.txt
              calculate_displacement 0_interface.dat $(ls -v *_interface.dat | tail -n 1)
              mkdir displacement_data
              mv *_displacement_data.dat displacement_data
              cd ../..
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
