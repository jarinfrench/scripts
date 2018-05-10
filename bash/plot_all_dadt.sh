read -p "Please enter the element: " element
if [[ "${element,,}" =~ ^(al|cu|fe)$ ]]; then
  element=${element,,}
  element=${element^}
else
  echo "$element is not recognized.  Please enter either Fe, Cu, or Al"
  exit 3
fi

#num_files=$(ls T*_[1-9]*.txt | awk -F '_' '{print $1"_"$2"_"$3"_"$4"_"$5"_"$6"_"$7}' | uniq | wc -l)
num_files=$(ls T*_[1-9]*.txt | awk -F '_' '{print $1"_"$2"_"$3"_"$4"_"$5}' | uniq | wc -l)
num=1
for i in $(ls T*_[1-9]*.txt | awk -F '_' '{print $1"_"$2"_"$3"_"$4"_"$5}' | uniq); do
  echo "Processing file $num of $num_files"
  temp=$(echo $i | awk -F '_' '{print $1}' | cut -c 2-)
  # rad_char=$(echo $i | awk -F '_' '{print $3}')
  # if [ "$rad_char" == "large" ]; then
  #   rad=100
  # elif [ "$rad_char" == "medium" ]; then
  #   rad=75
  # elif [ "$rad_char" == "small" ]; then
  #   rad=50
  # elif [ "$rad_char" == "tiny" ]; then
  #   rad=30
  # else
  #   echo "Error!"
  #   continue
  # fi
  rad=50
  if [ -a "${i}_average.txt" ]; then
    echo "Data set $i already processed."
    num=$(($num+1))
    continue
  fi
  ls -v ${i}_[1-9]* | xargs average_data.py
  gnuplot5 -e T=$temp -e "el='${element}'" -e r=$rad -e "basename='$i'" mobility.plt
  num=$(($num+1))
done
