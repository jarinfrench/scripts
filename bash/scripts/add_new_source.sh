#! /usr/bin/env bash

getDir() {
  echo "$1" | awk -F/ '{
    n=NF-'$2'+1;
    if(n < 1) exit;
    for (i=n; i <= NF; i++) {
      if (i == n) printf("%s", $i);
      else printf("/%s", $i);
    }
  }'
}

refdir="$1"
full_file_path=${2}
directory=$(dirname "${full_file_path}")
file=$(basename "${full_file_path}")

n=1
dir1=$(getDir ${directory} ${n}) # the parent directory of the new file
base=$(getDir ${refdir} 1) # The base directory where everything happens
if [[ "${dir1}" == "${base}" ]]; then
  echo "${file} - Taken from " >> "${refdir}/Sources.txt"
else
  check1=$(getDir ${directory} $((n + 1))) # the parent of the file's parent (the grandparent if you will)
  while [[ "${check1%%/*}" != "${base}" ]]; do
    ((n++))
    dir1=$(getDir ${directory} ${n}) # get the full path one more directory up
    tmp=$(getDir ${directory} $((n + 1))) # get the full path another directory up
    check1="${tmp%%/*}" # chop off all but the oldest parent
  done
  echo "${dir1}/${file} - Taken from " >> "${refdir}/Sources.txt"
fi
