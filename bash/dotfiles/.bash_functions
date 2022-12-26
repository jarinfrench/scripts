# compresses the specified file(s) into individual 7zip archives
compress() {
  if [ "$#" -eq 0 ]; then
    echo "Please enter the file(s) to compress"
    return 1
  else
    for i in "$@"; do
      if 7z a -m0=lzma -mx=9 "${i}".7z ${i}; then
        echo "Removing file ${i}"
        rm ${i}
      else
        echo "Error compressing file ${i}"
        return 2
      fi
    done
  fi
}

copy() {
  if [ -z "$1" ]; then
    echo "Copies a file directly to the clipboard from the terminal"
    echo "Usage: copy <filename>"
    return 1
  else
    cat "$1" | xclip -selection clipboard
  fi
}

extract() {
 if [ -z "$1" ]; then
    # display usage if no parameters given
    echo "Usage: extract <path/file_name>.<zip|rar|bz2|gz|tar|tbz2|tgz|Z|7z|xz|ex|tar.7z|tar.bz2|tar.gz|tar.xz>"
    echo "       extract <path/file_name_1.ext> [path/file_name_2.ext] [path/file_name_3.ext]"
    return 1
 else
    for n in $@
    do
      if [ -f "$n" ] ; then
          case "${n%,}" in
            *.tar.bz2|*.tar.gz|*.tar.xz|*.tbz2|*.tgz|*.txz|*.tar)
                         tar xvf "$n"       ;;
            *.tar.7z)    7z x -mmt=off ./"$n" -aos | tar xvf - -C . ;;
            *.lzma)      unlzma ./"$n"      ;;
            *.bz2)       bunzip2 ./"$n"     ;;
            *.rar)       unrar x -ad ./"$n" ;;
            *.gz)        gunzip ./"$n"      ;;
            *.zip)       unzip ./"$n"       ;;
            *.z)         uncompress ./"$n"  ;;
            *.7z|*.arj|*.cab|*.chm|*.deb|*.dmg|*.iso|*.lzh|*.msi|*.rpm|*.udf|*.wim|*.xar)
                         7z x -mmt=off ./"$n" -aos    ;; # the -aos skips extracting of existing files
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
  if [ -z "$1" ]; then
    num_lines=10
  else
    num_lines=$1
  fi
    history | awk '{CMD[$5]++;count++;}END { for (a in CMD)print CMD[a] " " CMD[a]/count*100 "% " a;}' | grep -v "./" | column -c3 -s " " -t | sort -nr | nl |  head -n ${num_lines}
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
    echo "Using prompt: <user>:<dirpath> (git branch) [time]$"
    PS1="\[\033[01;32m\]\u\[\033[00m\]:\[\033[01;34m\]$(shorten_path.py)\[\033[01;33m\]$(__git_ps1)\[\033[00m\] \[\033[01;31m\][\T]\[\033[00m\]$ "
  elif [[ $1 -eq "2" ]]; then
    # My prompt with the Conda environment indicated
    echo "Using prompt: (conda env) <user>:<dirpath> (git branch) [time]$)"
    PS1="\[\033[33m\]${CONDA_PROMPT_MODIFIER}\[\033[00m\]\[\033[01;32m\]\u\[\033[00m\]:\[\033[01;34m\]$(shorten_path.py)\[\033[01;33m\]$(__git_ps1)\[\033[00m\] \[\033[01;31m\][\T]\[\033[00m\]$ "
  elif [[ $1 -eq "3" ]]; then
    # A prompt with only the $ symbol
    echo "Using prompt: $ "
    PS1="\033[01;32m$\033[00m "
  elif [[ $1 -eq "4" ]]; then
    # A prompt with the current directory
    echo "Using prompt: <current directory>$ "
    PS1="\[\033[01;34m\]$(shorten_path.py)\033[01;32m $\033[00m "
  elif [[ $1 -eq "5" ]]; then
    # A full prompt with user@host:<path>
    echo "Using prompt: <user>@<host>:<current directory>$ "
    PS1="\033[01;32m\u@\H:\[\033[01;34m\]$(shorten_path.py)\033[01;32m $\033[00m "
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
  if [ -z "$1" ]; then
    echo "Usage: parent <file> <N>"
    echo "       N is the Nth parent from file (default 1)"
    return 1
  fi

  file=$1
  if [ -z "$2" ]; then
    num_parent=1
  else
    num_parent=$2
  fi


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

runtime() {
  if [ -z "$1" ]; then
    echo "Usage: runtime <program_name(s)>"
    echo ""
    return 1
  fi

  processes=${*:1}
  ps -acxo etime,command | grep -- "${processes// /\|}"
}

torange() {
  while read num; do
    if [[ -z ${first} ]]; then
      first=${num}
      last=${num}
      continue
    fi
    if [[ num -ne $((last + 1)) ]]; then
      if [[ first -eq last ]]; then
        echo ${first}
      else
        echo ${first}-${last}
      fi
      first=${num}
      last=${num}
    else
      : $((last++))
    fi
  done < $1
}

mklnk() {
  if [ -z "$2" ]; then
    echo "Makes a hard link of the file to the original"
    echo "Usage: ml <source_file> <destination_file>"
    echo ""
    return 1
  fi

  abs_path_src=$(get_abs_filename $1)
  abs_path_dest=$(get_abs_filename $2)

  ln ${abs_path_src} ${abs_path_dest}
  return 0
}

# An efficient way of searching the bashhub history - uses fzf (fuzzy finder)
hs() {
  if [ -z "$1" ]; then
    n=100
  else
    n=$1
  fi
  eval $(bh -n $n | fzf)
}

body() {
  # Print the header (first line of input), then run the specified command on the body
  # Use in a pipeline, e.g. ps | body grep somepattern
  IFS= read -r header
  printf '%s\n' "${header}"
  "$@"
}

newfile() {
  # returns the name of a unique file with the specified base name
  # Example: newfile tmp.txt --> output tmp(_#)?.txt, where _# only appears
  # if there is already a file name tmp.txt in the current directory. #
  # increases by one for each file already found (e.g. tmp_1.txt, tmp_2.txt, etc.)
  if [[ -z "$1" ]]; then
    echo "Enter the name(s) of the file(s) you would like to check"
    return 1
  else
    for file in $@; do
      python3 -c "from myModules import verify_new_file; print(verify_new_file('"${file}"'))"
    done
  fi
}

wsclean() {
  # for each file passed in, delete leading whitespace, and extra whitespace in it
  if [[ -z "$1" ]]; then
    echo "Enter the name(s) of the file(s) to clean up"
    return 1
  else
    for file in "$@"; do
      sed -i 's/^[ \t]*//g;s/ \+/ /g' "${file}" # remove leading whitespace (1st substitution), change consecutive spaces to 1 space (2nd substitution)
    done
  fi
}

trim() {
  # from https://awesomeopensource.com/project/dylanaraps/pure-bash-bible#trim-all-white-space-from-string-and-truncate-spaces
  # Removes leading/trailing/extra spaces from a string
  # Usage: trim "string" "pattern"
  # if pattern is not passed, the default is space

  if [[ -z "$2" ]]; then
    pattern=" "
    set -f
    set -- $*
    printf '%s\n' "$*"
    set +f
  else
    pattern="$2"
    printf '%s\n' "${1//${pattern}}"
  fi  
}

hex2rgb() {
  # from https://awesomeopensource.com/project/dylanaraps/pure-bash-bible#convert-a-hex-color-to-rgb
  # Usage: hex_to_rgb "#FFFFFF"
  #        hex_to_rgb "000000"
  : "${1/\#}"
  ((r=16#${_:0:2},g=16#${_:2:2},b=16#${_:4:2}))
  printf '%s\n' "$r $g $b"
}

rgb2hex() {
  # from https://awesomeopensource.com/project/dylanaraps/pure-bash-bible#convert-an-rgb-color-to-hex
  # Usage: rgb_to_hex "r" "g" "b"
  printf '#%02x%02x%02x\n' "$1" "$2" "$3"
}

source ~/.config/up/up.sh # see README for where to get this file.
source ~/projects/scripts/bash/optparse.bash

if [ -n "${SSH_CLIENT}" ] || [ -n "${SSH_TTY}" ]; then
  # The following are from Cascades .bash_functions
  # Various bash functions defined by the user
  cdls(){
    cd $1 && ls
  }

  scancel_range() {
    if [ "$#" -ne 2 ]; then
      echo "Please enter the beginning and end range for the jobs to delete"
      return 1
    else
      local beg=$1
      local end=$2
    fi

    job_list=`sq | awk -F " " 'NR>=2 {print $1}' | awk -F "." -v start=${beg} -v finish=${end} '{if ($1 >= start && $1 <= finish) print $1}'`
    n_jobs=`sq | awk -F " " 'NR>=2 {print $1}' | awk -F "." -v start=${beg} -v finish=${end} '{if ($1 >= start && $1 <= finish) print $1}' | wc -l`
    read -p "Are you sure you want to delete $n_jobs jobs? " verify
    case "$verify" in
      Y|y|[Yy][Ee][Ss] )  for i in ${job_list}; do scancel ${i}; done ;;
      N|n|[Nn][Oo] ) echo "Not deleting jobs."; return 0 ;;
      * ) echo "Invalid option" ;;
    esac
  }

  shist() {
    if [[ "$*" == *help* ]]; then
      echo "Show the history of jobs you've submitted. The first positional parameter is the number of jobs to show, and the second (optional) parameter can be one of either -c, -t, -p, -r, or -x to show jobs that have completed successfully, timed out, are pending, are running, or were cancelled respectively." | fold -sw 80
      return 0
    fi

    if [ -z "$1" ]; then
      local n=5
    else
      local n=$1
    fi
    local date_str=`date --date 'today -1 month' +%F`
    local format_str="JobID%10,JobName%15,Elapsed%13,State%12,End%20"
    if [ -z "$2" ]; then
      local cmd=DEFAULT
    else
      case "$2" in
        c|-c|--c|--completed|--complete|C|-C|--C) local cmd=COMPLETED ;;
        t|-t|--t|--timeout|T|-T|--T) local cmd=TIMEOUT ;;
        p|-p|--p|--pending|P|-P|--P) local cmd=PENDING ;;
        r|-r|--r|--running|R|-R|--R) local cmd=RUNNING ;;
        x|-x|--x|--cancelled|X|-X|--X) local cmd=CANCELLED ;;
        d|-d|--d|--default|D|-D|--D) local cmd=DEFAULT ;;
        *) echo "Unrecognized option, showing default"
           local cmd=DEFAULT ;;
      esac
    fi

    if [ "${cmd}" = "DEFAULT" ]; then
      sacct --starttime ${date_str} --format=${format_str} | awk 'NR==1 {print $0} NR==2 {print $0} {if ($2 != "extern" && $2 != "batch" && $2 != "orted" && $4 != "PENDING") a[++i%n]=$0} END {for (j=0;j<n;j++) print a[++i%n]}' n=${n}
    elif [ "${cmd}" = "COMPLETED" ]; then
      sacct --starttime ${date_str} --format=${format_str} | awk 'NR==1 {print $0} NR==2 {print $0} {if ($2 != "extern" && $2 != "batch" && $2 != "orted" && $4 == "COMPLETED") a[++i%n]=$0} END {for (j=0;j<n;j++) print a[++i%n]}' n=${n}
    elif [ "${cmd}" = "TIMEOUT" ]; then
      sacct --starttime ${date_str} --format=${format_str} | awk 'NR==1 {print $0} NR==2 {print $0} {if ($2 != "extern" && $2 != "batch" && $2 != "orted" && $4 == "TIMEOUT") a[++i%n]=$0} END {for (j=0;j<n;j++) print a[++i%n]}' n=${n}
    elif [ "${cmd}" = "PENDING" ]; then
      sacct --starttime ${date_str} --format=${format_str} | awk 'NR==1 {print $0} NR==2 {print $0} {if ($2 != "extern" && $2 != "batch" && $2 != "orted" && $4 == "PENDING") a[++i%n]=$0} END {for (j=0;j<n;j++) print a[++i%n]}' n=${n}
    elif [ "${cmd}" = "RUNNING" ]; then
      sacct --starttime ${date_str} --format=${format_str} | awk 'NR==1 {print $0} NR==2 {print $0} {if ($2 != "extern" && $2 != "batch" && $2 != "orted" && $4 == "RUNNING") a[++i%n]=$0} END {for (j=0;j<n;j++) print a[++i%n]}' n=${n}
    elif [ "${cmd}" = "CANCELLED" ]; then
      sacct --starttime ${date_str} --format=${format_str} | awk 'NR==1 {print $0} NR==2 {print $0} {if ($2 != "extern" && $2 != "batch" && $2 != "orted" && $4 == "CANCELLED") a[++i%n]=$0} END {for (j=0;j<n;j++) print a[++i%n]}' n=${n}
    else
      echo "Error: unknown command ${cmd}"
      return 10
    fi
  }

  cdjob() {
    if [ "$#" -ne 1 ]; then
      echo "Please enter the number of the job you wish to change to the directory of"
      return 1
    else
      local job_num=$1
    fi

    local dir_value=$(\squeue --user=jarinf -O workdir:200 -j ${job_num} 2>/dev/null | awk 'NR==2 {print $0}')
    if [ -z ${dir_value} ]; then
      local dir_value=$(sacct -j ${job_num} --format="WorkDir%200" | head -n 3 | tail -n 1 | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]] *$//')
    fi

    local RED='\033[0;31m'
    local NC='\033[0m'
    echo -e "${RED}${dir_value}${NC}"

    cd ${dir_value}
  }

  torange() {
    # Reads a file containing a list of numbers and converts it to a series of
    # ranges
    while read num; do
      if [[ -z ${first} ]]; then
        first=${num}
        last=${num}
        continue
      fi
      if [[ num -ne $((last + 1)) ]]; then
        if [[ first -eq last ]]; then
          echo ${first}
        else
          echo ${first}-${last}
        fi
        first=${num}
        last=${num}
      else
        : $((last++))
      fi
    done < $1
  }


################################################################################
# And the following are from Falcon
################################################################################
  # change to a directory and list the contents
  cdls() {
  cd "$1" && ls
  }

  # Find the working directory and cd to it of a job
  cdjob() {
  if [ "$#" -ne 1 ]; then
    echo "Please enter the job number you want to find."
    return 1
  fi
  # Handles the multi-line directories.
  new_dir=`\qstat -fx $1 | sed -n '/PBS_O_WORKDIR=/{:a;N;/,PBS_O_LANG/!ba;s/[[:space:]]//g;s/.*PBS_O_WORKDIR=\|,PBS_O_LANG.*//g;p}'`
  if [[ -z "$new_dir" ]]; then # handles the one-line directories
    new_dir=`\qstat -fx $1 | awk -F '=|,' '/PBS_O_WORKDIR=/{print $2}/PBS_O_LANG/{next}'`
  fi

  cd ${new_dir}
  pwd
  }

  up() {
    cd $(eval "printf '../'%.0s {1..$1}") && pwd
  }

  qalter_range() {
  if [ "$#" -ne 2 ]; then
    echo "Please enter the beginning and end range for the jobs to alter"
    return 1
  else
    local beg=$1
    local end=$2
  fi

  read -p "Please enter the options you want to change: " options
  job_list=`qstat | awk -F " " 'NR>=6 {print $1}' | awk -F "." -v start=${beg} -v finish=${end} '{if ($1 >= start && $1 <= finish) print $1}'`

  for i in ${job_list}; do
    qalter $i ${options}
  done
  }

  qdel_range() {
  if [ "$#" -ne 2 ]; then
    echo "Please enter the beginning and end range for the jobs to delete"
    return 1
  else
    local beg=$1
    local end=$2
  fi

  job_list=`qstat | awk -F " " 'NR>=6 {print $1}' | awk -F "." -v start=${beg} -v finish=${end} '{if ($1 >= start && $1 <= finish) print $1}'`
  n_jobs=`qstat | awk -F " " 'NR>=6 {print $1}' | awk -F "." -v start=${beg} -v finish=${end} '{if ($1 >= start && $1 <= finish) print $1}' | wc -l`
  read -p "Are you sure you want to delete $n_jobs jobs? " verify
  case "$verify" in
    Y|y|[Yy][Ee][Ss] )  for i in ${job_list}; do qdel ${i}; done ;;
    N|n|[Nn][Oo] ) echo "Not deleting jobs."; return 0 ;;
    * ) echo "Invalid option" ;;
  esac
  }

  compress() {
  if [ "$#" -eq 0 ]; then
    echo "Please enter the file(s) to compress"
    return 1
  else
    for i in "$@"; do
      if 7z a "${i}".7z ${i}; then
        echo "Removing file ${i}"
        rm ${i}
      else
        echo "Error compressing file ${i}"
        return 2
      fi
    done
  fi
  }

  decompress() {
  if [ "$#" -eq 0 ]; then
    echo "Please enter the file(s) to decompress"
    return 1
  else
    for i in "$@"; do
      7z x "${i}"
    done
  fi
  }
fi
