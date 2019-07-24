extract() {
 if [ -z "$1" ]; then
    # display usage if no parameters given
    echo "Usage: extract <path/file_name>.<zip|rar|bz2|gz|tar|tbz2|tgz|Z|7z|xz|ex|tar.bz2|tar.gz|tar.xz>"
    echo "       extract <path/file_name_1.ext> [path/file_name_2.ext] [path/file_name_3.ext]"
    return 1
 else
    for n in $@
    do
      if [ -f "$n" ] ; then
          case "${n%,}" in
            *.tar.bz2|*.tar.gz|*.tar.xz|*.tbz2|*.tgz|*.txz|*.tar)
                         tar xvf "$n"       ;;
            *.lzma)      unlzma ./"$n"      ;;
            *.bz2)       bunzip2 ./"$n"     ;;
            *.rar)       unrar x -ad ./"$n" ;;
            *.gz)        gunzip ./"$n"      ;;
            *.zip)       unzip ./"$n"       ;;
            *.z)         uncompress ./"$n"  ;;
            *.7z|*.arj|*.cab|*.chm|*.deb|*.dmg|*.iso|*.lzh|*.msi|*.rpm|*.udf|*.wim|*.xar)
                         7z x ./"$n"        ;;
            *.xz)        unxz ./"$n"        ;;
            *.exe)       cabextract ./"$n"  ;;
            *)
                         echo "extract: '$n' - unknown archive method"
                         return 1
                         ;;
          esac
      else
          echo "'$n' - file does not exist"
          return 1
      fi
    done
fi
}

hist-check() {
  if [ -z $1 ]; then
    num_lines=10
  else
    num_lines=$1
  fi
    history | awk '{CMD[$2]++;count++;}END { for (a in CMD)print CMD[a] " " CMD[a]/count*100 "% " a;}' | grep -v "./" | column -c3 -s " " -t | sort -nr | nl |  head -n ${num_lines}
}

mcd() {
  if [ -z $1 ]; then
    echo "Usage: mcd <new_directory_path>"
    return 1
  fi

  mkdir -p $1
  cd $1
}

dp() {
  if [[ $1 -eq "1" || $# -eq "0" ]]; then
    # My prompt
    PS1="\[\033[36m\]\u\[\033[m\]@\[\033[32m\]\h:\[\033[33;1m\]\w\[\033[m\]\$(__git_ps1 '(%s)')$ "
  elif [[ $1 -eq "2" ]]; then
    # A prompt with only the $ symbol
    PS1="\033[01;32m$\033[00m "
  elif [[ $1 -eq "3" ]]; then
    # A prompt with the current directory
    PS1="${debian_chroot:+($debian_chroot)}\w\033[01;32m$\033[00m "
  elif [[ $1 -eq "4" ]]; then
    # A full prompt with user@host:<path>
    PS1="\033[01;32m\u@\H:${debian_chroot:+($debian_chroot)}\w\033[01;32m$\033[00m "
  fi
  return;
}

get_abs_filename() {
  # $1 : relative filename
  filename=$1
  parentdir=$(dirname "${filename}")

  if [ -d "${filename}" ]; then
      echo "$(cd "${filename}" && pwd)"
  elif [ -d "${parentdir}" ]; then
    echo "$(cd "${parentdir}" && pwd)/$(basename "${filename}")"
  fi
}

parent() {
  if [ -z "$2" ]; then
    echo "Usage: parent <file> <N>"
    echo "       N is the Nth parent from file"
    return 1
  fi

  file=$1
  num_parent=$2

  if [ ${num_parent} -gt 100 ]; then
    echo "Are you sure you have more than 100 subdirectories? I'm assuming not."
    num_parent=100
  fi

  if [ ${num_parent} -eq 0 ]; then
    echo ${file}
    return 0
  fi

  full_path=$(get_abs_filename ${file})
  parent_dir=$(echo ${full_path} | sed -e 's;\/[^/]*$;;') # cut away /file.txt
  num=1
  while [ ${num} -lt ${num_parent} ]; do
  #for i in $(seq 1 $num_parent}); do
    parent_dir=$(echo ${parent_dir} | sed -e 's;\/[^/]*$;;') # cut away the next parent directory
    if [ -z "${parent_dir}" ]; then
      parent_dir="\\"
      break
    fi
    ((num++))
  done

  parent_dir=$(echo ${parent_dir} | sed -e 's;.*\/;;') # cut off everything before the parent dir
  echo "${parent_dir}"
  return 0
}

ProgressBarInit() {
  printf "\n"
  _last_progress=-1
}

# ProgressBar function - Input is currentState($1), totalState($2) and optionally the processName($3)
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
  if [ -z "$3" ]; then
    printf "\rProgress: [${_fill// /#}${_empty// /-}] ${_progress}%%"
  else
    printf "\r$3 Progress: [${_fill// /#}${_empty// /-}] ${_progress}%%"
  fi
}

# To be used in conjunction with alias lln (defined in .bash_aliases)
lf() {
 if [ "x${1}" == "x" ]; then
   n=1
 else
   n="${1}"
 fi
 ls -rt1 | tail -n ${n} | head -n 1
}

up() {
  # default parameter to 1 if non provided
  declare -i d=${@:-1}
  # ensure given parameter is non-negative. Print error and return if it is
  (( $d < 0 )) && (>&2 echo "up: Error: negative value provided") && return 1;
  # remove last d directories from pwd, append "/" in case result is empty
  cd "$(pwd | sed -E 's;(/[^/]*){0,'$d'}$;;')/";
}
