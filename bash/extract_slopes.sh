#! /bin/bash

if [[ -z "$1" ]]; then
  echo "Indicate 'pure' or 'impure' with this command"
  exit
fi
if [[ "$1" == "pure" ]]; then
  for i in $(fd dir_ -t d); do
    a=($(echo ${i} | awk -F'/' '{for (i=1;i<=NF;i++) {print $i}}'));
    ax=${a[0]}; 
    if [[ ${a[1]} == "sigma7" ]]; then 
      mis=38.20; 
    else 
      mis=${a[1]%degree}; 
    fi; 
    if [[ ${a[2]} == "Cooper" ]]; then 
      pot="CRG";
      pnum=17
    else 
      pot=${a[2]};
      pnum=1
    fi; 
    T=${a[3]#T}; 
    if [[ ${a[4]%_r} == "large" ]]; then 
      r=100; 
    elif [[ ${a[4]%_r} == "medium" ]]; then 
      r=75; 
    elif [[ ${a[4]%_r} == "small" ]]; then 
      r=50; 
    else 
      r=30; 
    fi; 
    n=$(echo ${a[${#a[@]}-1]} | awk -F'_' '{print $NF}'); 
    echo "${ax} ${pot} 0.0 Pure ${mis} ${T} ${r} ${n} $(tail -n1 ${i}/slope_calc.txt | awk '{print $4}') $(mobility_plots.py ${i}/area_data.txt -g | awk '{print $NF}' | tr -d '()')"; 
  done
elif [[ "$1" == "impure" ]]; then
  for i in $(fd dir_ -t d); do
    a=($(echo ${i} | awk -F'/' '{for (i=1;i<=NF;i++) {print $i}}'));
    ax=${a[0]}; 
    pot=${a[1]};
    c=$(echo ${a[2]} | cut -d'%' -f1)
    imp=$(echo ${a[2]} | cut -d'%' -f2)
    if [[ ${a[3]} == "sigma7" ]]; then 
      mis=38.20; 
    else 
      mis=${a[3]%degree}; 
    fi; 
    T=${a[4]#T}; 
    if [[ ${a[5]%_r} == "large" ]]; then 
      r=100; 
    elif [[ ${a[5]%_r} == "medium" ]]; then 
      r=75; 
    elif [[ ${a[5]%_r} == "small" ]]; then 
      r=50; 
    else 
      r=30; 
    fi; 
    n=$(echo ${a[${#a[@]}-1]} | awk -F'_' '{print $NF}'); 
    echo "${ax} ${pot} ${c} ${imp} ${mis} ${T} ${r} ${n} $(tail -n1 ${i}/slope_calc.txt | awk '{print $4}') $(mobility_plots.py ${i}/area_data.txt -g | awk '{print $NF}' | tr -d '()')"; 
  done
fi