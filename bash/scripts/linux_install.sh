#! /bin/bash

linkDotFile() {
  if [[ -z "$2" ]]; then
    dest="${HOME}/${1}"
  else
    dest="${HOME}/${2}"
  fi
  dateStr=$(date +%Y-%m-%d-%H%M) # current date to the minute

  if [[ -h "${HOME}/${1}" ]]; then
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
  ln -s "~/projects/scripts/bash/dotfiles/${1}" "${dest}"
}

shopt -s expand_aliases
box=$(uname -s)
RED='\033[0;31m'
if ! command -v curl &> /dev/null; then
  echo "curl is required for this install, but is not found. Installing curl..."
  if [[ "${box}" == "Linux" ]]; then
    sudo apt-get install -y curl
  else
    exit 0
  fi
fi

GREEN='\e[0;32m'
NC='\033[0m'
alias has='curl -sL https://git.io/_has | HAS_ALLOW_UNSAFE=y bash -s'

if [[ "${box}" != "Linux" ]]; then
  exit
fi

if ! has apt; then
  echo -e "${RED}apt required for installation${NC}"
  exit
fi

echo -e "${GREEN}Linux system detected, installing packages via apt${NC}"
# prepare atom for installation
wget -qO - https://packagecloud.io/AtomEditor/atom/gpgkey | sudo apt-key add -
sudo sh -c 'echo "deb [arch=amd64] https://packagecloud.io/AtomEditor/atom/any/ any main" > /etc/apt/sources.list.d/atom.list'
apt_packages=""
for i in mr task cloc htop parallel gawk atom autojump make gnuplot shellcheck qalc; do # add fdfind to this list if using Ubuntu 19.04 or later
  if ! has ${i} > /dev/null; then
    if [[ "${i}" == "mr" ]]; then
      apt_packages="${apt_packages} myrepos"
    elif [[ "${i}" == "task" ]]; then
      apt_packages="${apt_packages} taskwarrior"
      task_setup=y
    elif [[ "${i}" == "autojump" ]]; then
      link_autojump=y
    elif [[ "${i}" == "qalc" ]]; then
      apt_packages="${apt_packages} qalculate"
    # elif [[ "${i}" == "fdfind" ]]; then
    #   apt_packages="${apt_packages} fd-find"
    else
      apt_packages="${apt_packages} ${i}"
    fi
  fi
done

if ! has pip3 > /dev/null; then
  apt_packages="${apt_packages} python3-pip"
fi

if [[ ! -z "${apt_packages}" ]]; then
  echo -e "${GREEN}The following will be installed via apt: ${apt_packages}${NC}"
  sudo apt update -y
  sudo apt install -y ${apt_packages}
fi

if [[ "${task_setup}" == "y" ]]; then
  echo -e "${GREEN}Configuring taskwarrior. Type yes when prompted${NC}"
  task config taskd.server freecinc.com:53589
  task config taskd.key ~/.task/freecinc_b77224d4.key.pem
  task config taskd.certificate ~/.task/freecinc_b77224d4.cert.pem
  task config taskd.ca ~/.task/freecinc_b77224d4.ca.pem
  task config taskd.credentials -- 'FreeCinc/freecinc_b77224d4/cfdcc1e7-359c-4291-b7cc-34f91ad4d8b6'
fi

if [[ "${link_autojump}" == "y" ]]; then
  linkDotFile autojump.txt .local/share/autojump/autojump.txt
fi
