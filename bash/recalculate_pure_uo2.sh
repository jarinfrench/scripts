#! /bin/bash
extract () 
{ 
    if [ -z "$1" ]; then
        echo "Usage: extract <path/file_name>.<zip|rar|bz2|gz|tar|tbz2|tgz|Z|7z|xz|ex|tar.7z|tar.bz2|tar.gz|tar.xz>";
        echo "       extract <path/file_name_1.ext> [path/file_name_2.ext] [path/file_name_3.ext]";
        return 1;
    else
        for n in $@;
        do
            if [ -f "$n" ]; then
                case "${n%,}" in 
                    *.tar.bz2 | *.tar.gz | *.tar.xz | *.tbz2 | *.tgz | *.txz | *.tar)
                        tar xvf "$n"
                    ;;
                    *.tar.7z)
                        7z x -mmt=off ./"$n" -aos | tar xvf - -C .
                    ;;
                    *.lzma)
                        unlzma ./"$n"
                    ;;
                    *.bz2)
                        bunzip2 ./"$n"
                    ;;
                    *.rar)
                        unrar x -ad ./"$n"
                    ;;
                    *.gz)
                        gunzip ./"$n"
                    ;;
                    *.zip)
                        unzip ./"$n"
                    ;;
                    *.z)
                        uncompress ./"$n"
                    ;;
                    *.7z | *.arj | *.cab | *.chm | *.deb | *.dmg | *.iso | *.lzh | *.msi | *.rpm | *.udf | *.wim | *.xar)
                        7z x -mmt=off ./"$n" -aos
                    ;;
                    *.xz)
                        unxz ./"$n"
                    ;;
                    *.exe)
                        cabextract ./"$n"
                    ;;
                    *)
                        echo "extract: '$n' - unknown archive method";
                        return 1
                    ;;
                esac;
            else
                echo "'$n' - file does not exist";
                return 1;
            fi;
        done;
    fi
}

cd /media/jarinf/Research1/uo2/grain_growth/cylindrical
for i in $(fd dir -t d); do
  (
  cd ${i}
  if [[ -f snapshots.tar.xz ]]; then
    extract snapshots.tar.xz
  else
    if (($(ls *.dump.7z 2> /dev/null | wc -l) > 0)); then
      if (( ! $(ls *.dump 2> /dev/null | wc -l) == $(ls *.dump.7z 2> /dev/null | wc -l) )); then
        extract *.dump.7z
      fi
      tar -cJvf snapshots.tar.xz *.dump && rm *.dump.7z
    elif [[ -f snapshots.7z ]]; then
      if (( ! $(7z l snapshots.7z | tail -n1 | awk '{print $(NF-1)}') == $(ls *.dump 2> /dev/null | wc -l) )); then
        extract snapshots.7z
      fi
      tar -cJvf snapshots.tar.xz *.dump && rm snapshots.7z
    elif [[ -f snapshots.tar.bz2 ]]; then
      if (( ! $(tar -jtvf snapshots.tar.bz2 | wc -l) == $(ls *.dump 2> /dev/null | wc -l) )); then
        extract snapshots.tar.bz2
      fi
      tar -cJvf snapshots.tar.xz *.dump && rm snapshots.tar.bz2
    else
      echo "Cannot find dump files in ${i}"
      continue
    fi
  fi

  ax=$(echo ${i} | awk -F'/' '{print $1}')
  mis=$(echo ${i} | awk -F'/' '{print $2}')
  if [[ ${mis} == "sigma7" ]]; then
    mis=38.20
  else
    mis=${mis%degree}
  fi
  pot=$(echo ${i} | awk -F'/' '{print $3}')
  if [[ ${pot} == "Basak" ]]; then
    pnum=1
  elif [[ ${pot} == "Cooper" ]]; then
    pnum=17
  else
    echo "Unknown potential ${pot} in ${i}"
    continue
  fi
  n=2
  rcut=1.207
  a0=5.454
  s="fcc"
  echo "${mis} ${n} ${rcut} ${a0} ${s}" > find_grains_input2.txt
  case ${ax} in
    100)
      echo -e "1 0 0\n0 1 0\n0 0 1" >> find_grains_input2.txt
      ;;  
    110)
      echo -e "0 0 1\n1 -1 0\n1 1 0" >> find_grains_input2.txt
      ;;  
    111)
      echo -e "1 -1 0\n1 1 -2\n1 1 1" >> find_grains_input2.txt
      ;;  
    112)
      echo -e "1 1 1\n1 -1 0\n1 1 -2" >> find_grains_input2.txt
      ;;  
    *)  
      echo "ERROR: Unrecognized axis ${axis}"
      exit
      ;;  
  esac
  ls -v *.dump >> find_grains_input2.txt
  if [[ -f data.txt ]]; then
    mv data.txt data_original.txt
  fi
  find_grains_v2 find_grains_input2.txt -i 2 -e 0 -n 0.0
  if [[ -f "area_data.txt" ]]; then
    mv area_data.txt area_data_original.txt
  fi
  calculate_grain_area grain_area_input.txt -p ${pnum}
  mobility_plots.py area_data.txt --simple -F
  )
done

cd ${OLDPWD}