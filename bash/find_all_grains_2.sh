#! /bin/bash

read -p "Please enter the element: " element
if [[ "${element,,}" =~ ^(al|cu|fe)$ ]]; then
  element=${element,,}
  element=${element^}
else
  echo "$element is not recognized.  Please enter either Fe, Cu, or Al"
  exit 3
fi

if [[ "${element}" =~ ^Fe$ ]]; then
  potential=5
  lattice_param=2.85532
  cut100=0.933
  cut110=0.933
  cut111=0.933
  cutoff=1.6
elif [[ "${element}" =~ ^Al$ ]]; then
  potential=4
  lattice_param=4.032
  cut100=0.866
  cut110=0.866
  cut111=1.2
  cutoff=1.25
elif [[ "${element}" =~ ^Cu$ ]]; then
  potential=3
  lattice_param=3.615
  cut100=0.866
  cut110=0.866
  cut111=1.2
  cutoff=1.25
else
  potential=0
  read -p "Please enter the lattice parameter of the material: " lattice_param
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
  for j in $(ls -vd T*/); do
    cd ${j}
    for k in $(ls -d dir*/); do
      cd ${k}
      echo -e "\nProcessing directory ${i}/${j}${k}"
      if [ -d "interfaces" ]; then
        echo "Directory \"interfaces\" already found in ${i}/${j}${k}"
        t=$(echo $j | cut -c 2-)
        height=$(cat ${element}_minimized_* | head -n 7 | tail -n 1 | awk '{print $2}')
        calculate_grain_area data.txt $t $height ${lattice_param} $potential
        calculate_force_and_velocity.py $t ${height} ${lattice_param} -p $potential
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
          t=$(echo $j | cut -c 2-)
          height=$(cat ${element}_minimized_* | head -n 7 | tail -n 1 | awk '{print $2}')
          calculate_grain_area data.txt $t ${height} ${lattice_param} $potential
          calculate_force_and_velocity.py $t ${height} ${lattice_param} -p $potential
        else
          echo -e "\033[0;32m\tArea data file already found in ${i}/${j}${k}\033[0m"
        fi
      fi
      cd ..
    done
    cd ..
  done
  cd ..
done
