#! /bin/bash

if ! command -v curl &> /dev/null; then
  echo "curl is required for this install, but is not found. Installing curl..."
  box=$(uname -s)
  if [[ "${box}" == "Linux" ]]; then
    sudo apt-get install curl
  elif [[ "${box}" == "Darwin" ]]; then
    if command -v brew &> /dev/null; then
      brew install curl
    else
      echo "Homebrew not found on ${box} system, unable to install curl. Exiting..."
      exit 1
    fi
  else
    echo "No installation procedures for ${box}, exiting..."
    exit 0
  fi
fi

bash -c ./install_dotfiles.sh
. ~/.bashrc
bash -c ./brew_install.sh
bash -c ./linux_install.sh
bash -c ./pip_install.sh
bash -c ./source_install.sh
bash -c ./npm_install.sh
bash -c ./snap_install.sh
bash -c ./git_install.sh
bash -c ./extras_install.sh
bash -c ./python_scripts_install.sh
