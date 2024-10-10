#! /usr/bin/env bash

dir2check="/media/jarinf/Research1/Images/Outside_sources"
image_regex='^(jpe?g|png|[g|t]iff?|bmp|eps)$'

# Note that the directory read in is the absolute directory
inotifywait -qmre create -e moved_to --format '%w %f %e' "${dir2check}" | while read directory file action; do
  if [[ "${action}" == "MOVED_TO,ISDIR" ]]; then
    # echo "Directory moved here"
    find "${directory}${file}" -type f -regextype posix-extended -regex '^.*\.(jpe?g|png|[g|t]iff?|bmp|eps)$' -exec /home/jarinf/projects/scripts/bin/add_new_source.sh "${dir2check}" {} \;
  elif [[ ! ${file##*.} =~ ${image_regex} ]]; then
    # echo "File ${file} (extension ${file##*.}) ignored"
    continue # ignore all files that are not images.
  elif [[ "${action}" == "CREATE,ISDIR" ]]; then
    # New directory created - do nothing
    continue
  elif [[ "${action}" == "CREATE" ]]; then
    # echo "Created a new file"
    /home/jarinf/projects/scripts/bin/add_new_source.sh "${dir2check}" "${directory}${file}"
  elif [[ "${action}" == "MOVED_TO" ]]; then
    # echo "File was moved here"
    /home/jarinf/projects/scripts/bin/add_new_source.sh "${dir2check}" "${directory}${file}" # TODO; change this so all I do is move the relevant entry (can I do that - requires moved_from in inotifywait?)
  fi
  gedit "${dir2check}/Sources.txt"
  sort "${dir2check}/Sources.txt" -o "${dir2check}/Sources.txt"
done
