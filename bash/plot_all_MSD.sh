cd /media/jarinf/Seagate\ Expansion\ Drive/UO2/grain_growth/MSD_data
for i in $(ls T*_with_xy.dat); do
  temp=$(echo $i | awk -F '_' '{print $1}' | cut -c 2-)
  rad_char=$(echo $i | awk -F '_' '{print $3}')
  if [ "${rad_char}" == "large" ]; then
    rad=100
  elif [ "${rad_char}" == "medium" ]; then
    rad=75
  elif [ "${rad_char}" == "small" ]; then
    rad=50
  elif [ "${rad_char}" == "tiny" ]; then
    rad=30
  else
    echo "Error!"
    exit 2
  fi
  gnuplot -e T=$temp -e "el='UO2'" -e r=$rad -e "filename='$i'" MSD.plt
done
