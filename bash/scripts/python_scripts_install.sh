#! /bin/bash

pythonscriptspath="${HOME}/projects/scripts/python"
bindir="${HOME}/.local/bin"

linkPythonFile() {
  # $1 is the original full-path file name
  # $2, if given, is the desired link name
  if [[ -z "$2" ]]; then
    dest="${bindir}/${1}"
  else
    dest="${bindir}/${2}"
  fi
  dateStr=$(date +%Y-%m-%d-%H%M) # current date to the minute

  if [[ -h "${dest}" ]]; then
    # existing symlink
    echo "Removing existing symlink: ${dest}"
    rm "${dest}"
  elif [[ -f "${dest}" ]]; then
    # existing file
    echo "Backing up existing file: ${dest}"
    mv "${dest}"{,.${dateStr}}
  elif [[ -d "${dest}" ]]; then
    # existing directory
    echo "Backing up existing directory: ${dest}"
    mv "${dest}"{,.${dateStr}}
  fi

  echo "Creating new symlink: ${dest}"
  ln -s "${1}" "${dest}"
}

while IFS= read -r -d '' file
do
  name=$(basename "${file}")
  linkPythonFile "${file}" "${name}"
done < <(find "${pythonscriptspath}" -executable -type f -name "*.py" -print0)
