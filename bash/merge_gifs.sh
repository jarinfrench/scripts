#! /bin/bash

trap ctrl_c INT

ctrl_c() {
  exit 3
}

ProgressBarInit() {
  printf "\n"
  _last_progress=-1
}

# ProgressBar function - Input is currentState($1), totalState($2) and optional processName($3)
ProgressBar() {
  _progress=$((${1} * 100 / ${2}))

  if [ ${_progress} -eq "${_last_progress}" ]; then
    return
  fi

  _last_progress=${_progress}

  _done=$((_progress * 4 / 10))
  _left=$((40 - _done))

  # Build progressbar string lengths
  _fill=$(printf "%${_done}s")
  _empty=$(printf "%${_left}s")

  # Build progressbar strings and print the progressbar line
  if [ -z "$3" ]; then
    printf "\rProgress: [${_fill// /#}${_empty// /-}] ${_progress}%%"
  else
    printf "\r$3 Progress: [${_fill// /#}${_empty// /-}] ${_progress}%%"
  fi
}

#
# Take multiple *.gif animations and stack them within a 3x3 grid.  Assumes that
# the number of frames, the animation time, etc. are all equal, and that only
# the frames are different between the gifs.
# Arguments: base names of the original gifs (e.g. w/o the .gif extension)
# See http://www.imagemagick.org/Usage/anim_mods/#append

# Assume a resolution of 1920 x 1000
# max width of each frame should be 600 - total of 1800 across three images
# max height should then be 300 - total of 900 across three images

num_args=$#

if [ "${num_args}" -lt 2 ]; then
  echo "At least two gifs must be specified."
  exit 1
fi

if [ "${num_args}" -gt 9 ]; then
  echo "Too many gifs - unable to generate merged gif."
  exit 2
fi

case ${num_args} in
1)
  echo "At least two gifs must be specified."
  exit 1
  ;;
2)
  append_plus=("1 2")
  append_minus=()
  ;;

3)
  append_plus=("1 2" "comb12 3") # two separate commands - one to combine 1 and 2 together, then another command to combine that combined image with 3
  append_minus=()
  ;;

4)
  append_plus=("1 2" "3 4")      # two separate commands
  append_minus=("comb12 comb34") # combined 1 and 2, combined 3 and 4
  ;;

5)
  append_plus=("1 2" "comb12 3" "4 5")
  append_minus=("comb123" "comb45")
  ;;

6)
  append_plus=("1 2" "comb12 3" "4 5" "comb45 6")
  append_minus=("comb123 comb456")
  ;;

7)
  append_plus=("1 2" "comb12 3" "4 5" "comb45 6")
  append_minus=("comb123 comb456" "comb123456 7")
  ;;

8)
  append_plus=("1 2" "comb12 3" "4 5" "comb45 6" "7 8")
  append_minus=("comb123 comb456" "comb123456 comb78")
  ;;

9)
  append_plus=("1 2" "comb12 3" "4 5" "comb45 6" "7 8" "comb78 9")
  append_minus=("comb123 comb456" "comb123456 comb789")
  ;;

*)
  echo "Too many gifs - unable to generate merged gif."
  exit 2
  ;;
esac

# get each gif's individual frames with animation information
mkdir .tmp
cd .tmp || {
  echo "Failed to change to '.tmp'"
  exit
}
ProgressBarInit
iter=1
total=${num_args}
num_frames=()
ProgressBar 0 ${total} "Frame Splitting"
for _gif in "$@"; do
  # explode the gif into it's individual frames
  gif2anim -c ../"${_gif}".gif

  # resize and title each one.
  num_frames[${#num_frames[@]}]=$(ls ${_gif}_[0-9][0-9][0-9].gif | wc -l)
  for frame in "${_gif}"_[0-9][0-9][0-9].gif; do
    basename="${frame%.*}"
    new_name="${basename}_small.gif"
    convert "${frame}" -resize 900x900 "${new_name}"
    # This can be fiddled with to get better looking labels/titles, but it works for now.
    #convert "${new_name}" label:"${_gif}" +swap -gravity Center -append "${new_name}"
  done

  ProgressBar ${iter} ${total} "Frame Splitting"
  ((iter++))
done
if [[ "${#num_frames[@]}" -gt 0 ]] && [[ $(printf "%s\000" "${num_frames[@]}" | LC_ALL=C sort -z -u | grep -z -c .) -eq 1 ]]; then
  n_frames=${num_frames[0]}
else
  echo ""
  echo "Warning! Number of frames not consistent between gifs! Using minimum number of frames."
  min=${num_frames[0]}
  idx=0
  for val in "${num_frames[@]}"; do
    echo "Frames in gif ${idx}: ${val}"
    ((idx++))
    ((val < min)) && min=${val}
  done
  n_frames=${min}
fi

ProgressBarInit
iter=1
# take the number of strings in append_plus, added to the number of strings in append_minus
# and multiply by the number of frames.
total=$(((${#append_plus[@]} + ${#append_minus[@]}) * n_frames))
ProgressBar 0 ${total} "Frame Stitching (+)"
for pair in "${append_plus[@]}"; do
  str_arr=(${pair})
  index1=${str_arr[0]}
  index2=${str_arr[1]}
  re='^[1-9]$'
  if ! [[ ${index1} =~ ${re} ]]; then      # not a number, i.e. comb12
    new_name="${index1}${index2//[!1-9]/}" # also handles the case where index2 is not a number (i.e. comb34) - would end up with comb1234
  else
    new_name="comb${index1}${index2}" # assumes that index2 will always be a number if index1 is a number
  fi

  for frame in $(seq -f '%03g' 1 "${n_frames}"); do
    if ! [[ ${index1} =~ ${re} ]]; then   # not a number, i.e. comb12
      if ! [[ ${index2} =~ ${re} ]]; then # not a number, i.e. comb34
        convert "${index1}_${frame}.gif" "${index2}_${frame}.gif" +append "${new_name}_${frame}.gif"
      else
        convert "${index1}_${frame}.gif" "${!index2}_${frame}_small.gif" +append "${new_name}_${frame}.gif"
      fi
    else
      convert "${!index1}_${frame}_small.gif" "${!index2}_${frame}_small.gif" +append "${new_name}_${frame}.gif"
    fi

    ProgressBar ${iter} ${total} "Frame Stitching (+)"
    ((iter++))
  done
done

for pair in "${append_minus[@]}"; do
  str_arr=(${pair})
  index1=${str_arr[0]}
  index2=${str_arr[1]}
  re='^[1-9]$'
  # Note that for the -append option, index1 will _always_ be non-numeric
  new_name="${index1}${index2//[!1-9]/}"

  for frame in $(seq -f '%03g' 1 "${n_frames}"); do
    if ! [[ ${index2} =~ ${re} ]]; then # not a number, i.e. comb34
      convert "${index1}_${frame}.gif" "${index2}_${frame}.gif" -append "${new_name}_${frame}.gif"
    else
      convert "${index1}_${frame}.gif" "${!index2}_${frame}_small.gif" -append "${new_name}_${frame}.gif"
    fi
    ProgressBar ${iter} ${total} "Frame Stitching (-)"
    ((iter++))
  done
done

anim2gif -c -b "${new_name}" $1.anim
mv "${new_name}_anim.gif" ../combined.gif
echo ""
# cleanup
cd ../ || {
  echo "Failed to change to parent directory"
  exit 1
}
rm -rf .tmp
