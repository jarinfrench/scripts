trap ctrl_c INT

function ctrl_c() {
  # Force exiting of the script, rather than just the current executable
  exit 1
}

read -p "Please enter the element: " element
if [[ "${element,,}" =~ ^(al|cu|fe)$ ]]; then
  element=${element,,}
  element=${element^}
  num_files=$(ls T*_[1-9]*.txt | awk -F '_' '{print $1"_"$2"_"$3"_"$4"_"$5}' | uniq | wc -l)
  files=$(ls T*_[1-9]*.txt | awk -F '_' '{print $1"_"$2"_"$3"_"$4"_"$5}' | uniq)
elif [[ "${element,,}" =~ ^uo2$ ]]; then
  element=${element^^}
  num_files=$(ls T*_[1-9].txt | awk -F '_' '{print $1"_"$2"_"$3"_"$4"_"$5"_"$6"_"$7"_"$8}' | uniq | wc -l)
  files=$(ls T*_[1-9].txt | awk -F '_' '{print $1"_"$2"_"$3"_"$4"_"$5"_"$6"_"$7"_"$8}' | uniq)
elif [[ "${element,,}" =~ ^u$ ]]; then
  element=${element^^}
  num_files=$(ls T*_[1-9]*.txt | awk -F '_' '{print $1"_"$2"_"$3"_"$4"_"$5"_"$6"_"$7"_"$8"_"$9"_"$10}' | uniq | wc -l)
  files=$(ls T*_[1-9]*.txt | awk -F '_' '{print $1"_"$2"_"$3"_"$4"_"$5"_"$6"_"$7"_"$8"_"$9"_"$10}' | uniq)
else
  echo "$element is not recognized.  Please enter either UO2, Fe, Cu, U, or Al"
  exit 3
fi

num=1
for i in $files; do
  echo ${i}
  echo ""
  echo "Processing file $num of $num_files"
  temp=$(echo $i | awk -F '_' '{print $1}' | cut -c 2-)
  if [[ "${element}" =~ ^(UO2)$ ]]; then
    rad_char=$(echo $i | awk -F '_' '{print $4}')
    if [ "$rad_char" == "large" ]; then
      rad=100
    elif [ "$rad_char" == "medium" ]; then
      rad=75
    elif [ "$rad_char" == "small" ]; then
      rad=50
    elif [ "$rad_char" == "tiny" ]; then
      rad=30
    else
      echo "Error! Radius value not recognized."
      continue
    fi
  else
    rad=100
  fi
  if [ -a "${i}_average.txt" ]; then
    echo "Data set $i already processed."
    num=$(($num+1))
    continue
  fi
  xenon_amount=$(echo ${i} | awk -F '_' '{print $4}')
  ls -v ${i}_[1-9]* | xargs average_data.py
  gnuplot -e T=$temp -e "el='${element}+${xenon_amount}Xe'" -e r=$rad -e "basename='$i'" mobility.plt
  num=$(($num+1))
done
