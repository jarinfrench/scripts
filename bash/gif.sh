#! /bin/bash

ProgressBarInit() {
  printf "\n"
  _last_progress=-1
}

# ProgressBar function - Input is currentState($1) and totalState($2)
ProgressBar() {
  # Process data
  let _progress=(${1}*100/${2}*100)/100

  if [ ${_progress} -eq ${_last_progress} ]; then
    return
  fi

  _last_progress=${_progress}

  let _done=(${_progress}*4)/10
  let _left=40-${_done}

  # Build progressbar string lengths
  _fill=$(printf "%${_done}s")
  _empty=$(printf "%${_left}s")

  # Build progressbar strings and print the progressbar line
  printf "\rProgress: [${_fill// /#}${_empty// /-}] ${_progress}%%"
}

#
# Take multiple *.png sequences outputted by paraview and generate a
# single side by side gif animation
# Arguments: base names of the image series
#

# output file name
outname=`echo $@ | sed 's/ /_/g'`.gif

# montage the series together
mkdir .tmp
ProgressBarInit
iter=1
tot=`ls $1.*.png | wc -l`
for file in $1.*.png
do
  # get number
  num=`echo $file | cut -d. -f2`
  arg=''
  for base in $@
  do
    arg=$arg' '${base}.${num}.png
  done
  convert +append $arg .tmp/${num}.gif
  ProgressBar ${iter} ${tot}
  ((iter++))
done

# animate
for anim in `ls | cut -d. -f1|sort -u`
do
  gifsicle -d20 .tmp/*.gif --colors 256 -l > $outname 2>> errors.log
done

echo ""
# cleanup
rm -rf .tmp
