#! /bin/bash

shopt -s expand_aliases
dotfilesdir="${HOME}/projects/scripts/bash/dotfiles"
box=$(uname -s)

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
  ln -s "${dotfilesdir}/${1}" "${dest}"
}

linkDotFile .atom
linkDotFile .ssh
linkDotFile .task
linkDotFile .vim
linkDotFile .alias_completion.sh
linkDotFile .bash_aliases
linkDotFile .bash_command_timer.sh
linkDotFile .bash_functions
linkDotFile .bash_variables
linkDotFile .dircolors .dir_colors
linkDotFile .gdbinit
linkDotFile .gitconfig
linkDotFile .gnuplot
linkDotFile .inputrc
linkDotFile .jrnl_config
linkDotFile .mrconfig
linkDotFile .pdbrc
linkDotFile .taskrc
linkDotFile .vimrc


if [[ "${box}" == "Darwin" ]]; then
  linkDotFile .bash_profile
  linkDotFile .bashrc_mac .bashrc
elif [[ "${box}" == "Linux" ]]; then
  linkDotFile .bashrc_linux .bashrc
fi

if [[ ! -d "~/.Mathematica/Kernel" ]]; then
  mkdir -p ~/.Mathematica/Kernel
fi
linkDotFile init.m .Mathematica/Kernel/init.m
