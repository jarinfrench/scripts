read -p "Please enter the element: " element
if [[ "${element,,}" =~ ^(al|cu|fe)$ ]]; then
  element=${element,,}
  element=${element^}
  cd /media/jarinf/Seagate\ Expansion\ Drive/${element}/grain_growth/force_velocity_data
  files=$(ls T*_[1-9]*.txt | awk -F '_' '{print $1"_"$2"_"$3"_"$4"_"$5"_"$6}' | uniq)
  num_files=$(ls T*_[1-9]*.txt | awk -F '_' '{print $1"_"$2"_"$3"_"$4"_"$5"_"$6}' | uniq | wc -l)
elif [[ "${element,,}" =~ ^uo2$ ]]; then
  element=${element,,}
  element=${element^^}
  cd /media/jarinf/Seagate\ Expansion\ Drive/${element}/grain_growth/force_velocity_data
  files=$(ls T*_[1-9]*.txt | awk -F '_' '{print $1"_"$2"_"$3"_"$4"_"$5"_"$6"_"$7"_"$8}' | uniq)
  num_files=$(ls T*_[1-9]*.txt | awk -F '_' '{print $1"_"$2"_"$3"_"$4"_"$5"_"$6"_"$7"_"$8}' | uniq | wc -l)
else
  echo "$element is not recognized.  Please enter either UO2, Fe, Cu, or Al"
  exit 3
fi

num=1
for i in ${files}; do
  echo "Processing file $num of $num_files"
  temp=$(echo $i | awk -F '_' '{print $1}' | cut -c 2-)
  if [[ "${element}" ~= ^UO2$ ]]; then
    rad_char=$(echo $i | awk -F '_' '{print $3}')
    if [ "$rad_char" == "large" ]; then
      rad=100
    elif [ "$rad_char" == "medium" ]; then
      rad=75
    elif [ "$rad_char" == "small" ]; then
      rad=50
    elif [ "$rad_char" == "tiny" ]; then
      rad=30
    else
      echo "Error!"
      continue
    fi
  elif [[ "${element,,}" ~= ^(al|cu|fe)$ ]]; then
    rad=$(echo $i | awk -F '_' '{print $3}' | cut -c 2-) # At the moment, this is only 50, but this should change, and I may change the file name formatting of the Al, Cu, and Fe files.
  else
    echo "Something weird happened.  What element are we using again?"
    exit 4
  fi

  if [ -a "${i}_average.txt" ]; then
    echo "Data set $i already processed."
    num=$(($num+1))
    continue
  fi

  ls -v ${i}_[1-9]* | xargs average_data.py
  gnuplot -e T=$temp -e "el='${element}'" -e r=$rad -e "filename='$i'" plot_vel_vs_force.plt
