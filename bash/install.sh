#! /bin/bash

dotfiles="bash_aliases bash_functions bash_variables gdbinit gitconfig gnuplot inputrc jrnl_config mrconfig pdbrc taskrc vimrc"

function ctrl_c() {
  echo -e "Cleaning up..."
  for i in ${dotfiles} bashrc bash_profile; do
    if [[ -f ${HOME}/.${i}.bak ]]; then
      mv ${HOME}/.${i}.bak ${HOME}/.${i}
    fi
  done
  exit 1
}

alias has='curl -sL https://git.io/_has | HAS_ALLOW_UNSAFE=y bash -s'

GREEN='\033[0;32m'
RED='\033[0;31m'
NC='\033[0m'

trap ctrl_c INT

INSTALL_DIR="${HOME}/projects/scripts/bash"
LINK_DIR="${HOME}"
BIN_DIR="${HOME}/projects/scripts/bin"

cd ${INSTALL_DIR}

# -b uses a simple backup for the destination file (if it exists), and -S changes the suffix from ~ to .bak (in this case)
for i in ${dotfiles}; do
  if [ -f ${LINK_DIR}/.${i}.bak ]; then
    echo "${HOME}/.${i} already has a backup file (${HOME}/.${i}.bak) - unable to link file."
  else
    ln -sbS .bak ${INSTALL_DIR}/.${i} ${LINK_DIR}/.${i}
  fi
done

# ssh config
if [ -f ${LINK_DIR}/.ssh/config.bak ]; then
  echo "${RED}ssh config file already has a backup (${LINK_DIR}/.ssh/config.bak) - unable to link config file${NC}"
else
  ln -sbS .bak ${INSTALL_DIR}/ssh_config ${LINK_DIR}/.ssh/config
fi

# atom configuration file
if [ -f ${LINK_DIR}/.atom/config.cson.bak ]; then
  echo "${RED}Atom configuration file (config.cson) already has a backup (${LINK_DIR}/.atom/config.cson.bak)- unable to link config file${NC}"
else
  ln -sbS .bak ${INSTALL_DIR}/atom_config ${LINK_DIR}/.atom/config.cson
fi

ln -sf ${INSTALL_DIR}/.alias_completion.sh ${LINK_DIR}/.alias_completion.sh

# Store symlinks to python programs in bin
for i in $(find ~/projects/scripts/python/ -type f -executable); do
  file=$(basename ${i})
  ln -sf ${i} ${BIN_DIR}/${file}
done

# soft link the executables in the bash directory (excluding install.sh) to the bin directory
for i in $(find ~/projects/scripts/bash/ -maxdepth 1 -type f ! -iname "install.sh" -executable); do
  file=$(basename ${i})
  ln -sf ${i} ${BIN_DIR}/${file%.*} # links the executable bash file to the bin dir with the same name, minus the extension.
done

ln -sbS .bak ${INSTALL_DIR}/.vim ${LINK_DIR}/.vim

UNAME=$(uname -s)
if [[ "${UNAME}" == "Darwin" ]]; then
  echo -e "${GREEN}If you would like your current .bash_profile overwritten, execute the following command:"
  echo -e "ln -sbS .bak ${INSTALL_DIR}/.bash_profile ~/.bash_profile"
  echo -e "${NC}"

  echo -e "${GREEN}If you would like your current .bashrc overwritten, execute the following command:"
  echo -e "ln -sbS .bak ${INSTALL_DIR}/.bashrc_mac ~/.bashrc"
  echo -e "${NC}"

  # required programs:
  has snap curl brew wget npm git
  if [ $? -gt 0 ]; then
    echo -e "${RED}Make sure you have the above software installed!${NC}"
    exit 2
  fi
elif [[ "${UNAME}" == "Linux" ]]; then
  echo -e "${GREEN}If you would like your current .bashrc overwritten, execute the following command:"
  echo -e "ln -sbS .bak ${INSTALL_DIR}/.bashrc_linux ~/.bashrc"
  echo -e "${NC}"

  has snap curl wget npm git pip3
  if [ $? -gt 0 ]; then
    echo -e "${RED}Make sure you have the above software installed!${NC}"
    exit 2
  fi
fi

if [[ "${UNAME}" == "Darwin" ]]; then
  echo -e "Mac system detected, installing packages using Homebrew."
  # install using Darwin methods
  # Note that these methods have been tested on Sierra, and found to work!
  brew_packages=""
  for i in fd hr mr jrnl task taskd tasksh rng cloc fzf; do
    has ${i} > /dev/null
    if [ $? -gt 0 ]; then
      if [ "${i}" == "rng" ]; then
        TAP_BREW=y
      fi
      brew_packages="${brew_packages} ${i}"
    fi
  done

  echo -e "Installing the following packages from Homebrew: ${brew_packages}"
  has rng > /dev/null
  if [ "${TAP_BREW}" == "y" ]; then
    brew tap nickolasburr/pfa
  fi
  brew install ${brew_packages}

  source_packages="optparse.bash up"
  for i in bd has ansi bashhub; do
    has ${i} > /dev/null
    if [ $? -gt 0 ]; then
      source_packages="${source_packages} ${i}"
    fi
  done
  echo -e "Installing the following packages from source: ${source_packages}"

elif [[ "${UNAME}" == "Linux" ]]; then
  echo -e "${GREEN}Linux system detected, installing packages via apt."
  # install using UBUNTU methods (may need to change this later, but it should work for now)
  apt_packages=""
  for i in mr task cloc htop parallel; do
    has ${i} > /dev/null
    if [ $? -gt 0 ]; then
      if [ "${i}" == "mr" ]; then
        apt_packages="${apt_packages} myrepos"
      elif [ "${i}" == "task" ]; then
        apt_packages="${apt_packages} taskwarrior"
        TASK_SETUP=y
      else
        apt_packages="${apt_packages} ${i}"
      fi
    fi
  done
  echo -e "Installing the following packages from repositories: ${apt_packages}${NC}"
  sudo apt install ${apt_packages}

  # set up task files
  if [ "${TASK_SETUP}" == "y" ]; then
    echo -e "  ${GREEN}Configuring taskwarrior... please type yes when prompted${NC}"
    ln -sf ${INSTALL_DIR}/.task ${HOME}/.task
    task config taskd.server freecinc.com:53589
    task config taskd.key ~/.task/freecinc_b77224d4.key.pem
    task config taskd.certificate ~/.task/freecinc_b77224d4.cert.pem
    task config taskd.ca ~/.task/freecinc_b77224d4.ca.pem
    task config taskd.credentials -- 'FreeCinc/freecinc_b77224d4/cfdcc1e7-359c-4291-b7cc-34f91ad4d8b6'
  fi

  has jrnl > /dev/null
  if [ $? -gt 0 ]; then
    echo -e "${GREEN}Installing the following package using pip: jrnl${NC}"
    pip3 install --user jrnl
  fi

  source_packages="optparse.bash up"
  for i in bd fd has hr htop ansi bashhub fzf; do
    has ${i} > /dev/null
    if [ $? -gt 0 ]; then
      source_packages="${source_packages} ${i}"
    fi
  done
  echo -e "${GREEN}Installing the following packages from source: ${source_packages}"

  has fd > /dev/null
  if [ $? -gt 0 ]; then
    echo -e "  Installing fd${NC}"
    curl https://github.com/sharkdp/fd/releases/download/v7.3.0/fd-musl_7.3.0_amd64.deb > fd-musl_7.3.0_amd64.deb
    sudo dpkg -i fd-musl_7.3.0_amd64.deb
    rm fd-musl_7.3.0_amd64.deb
  fi

  has hr > /dev/null
  if [ $? -gt 0 ]; then
    echo -e "${GREEN}  Installing hr${NC}"
    sudo curl https://raw.githubusercontent.com/LuRsT/hr/master/hr > hr
    sudo chmod +x hr
    sudo mv hr /usr/local/bin/
  fi

  has rng > /dev/null
  if [ $? -gt 0 ]; then
    if [ ! -d ${HOME}/projects/scripts/bash/rng ]; then
      echo -e "${GREEN}  Installing rng${NC}"
      git clone https://github.com/nickolasburr/rng.git
      cd rng
      make
      sudo make install
      cd ../
    fi
  fi

  has fzf > /dev/null
  if [ $? -gt 0 ]; then
    if [ ! -d ${HOME}/.fzf ]; then
      echo -e "${GREEN}  Installing fzf${NC}"
      git clone --depth 1 https://github.com/junegunn/fzf.git ${HOME}/.fzf
      ${HOME}/.fzf/install
    fi
  fi
fi

has bashhub > /dev/null
if [ $? -gt 0 ]; then
  echo -e "${GREEN}  Installing bashhub${NC}"
  curl -OL https://bashhub.com/setup && bash setup # This requires user interaction
  rm setup
fi

. .bash_aliases # required to check if bd is aliased
has bd > /dev/null
if [ $? -gt 0 ]; then
  echo -e "${GREEN}  Installing bd${NC}"
  sudo wget --no-check-certificate -O /usr/local/bin/bd https://raw.github.com/vigneshwaranr/bd/master/bd
  sudo chmod +rx /usr/local/bin/bd
  alias bd 2>/dev/null >/dev/null || (echo -e 'alias bd=". bd -si"' >> .bash_aliases && source .bash_aliases) # checks if the alias bd exists, and if not, adds the alias to the alias list.
fi

echo -e "${GREEN}  Installing optparse.bash"
echo -e "    Note that optparse requires the GNU version of sed (for Mac - install by brew install gnu-sed (--with-default-names, if you don't want to alias sed))${NC}"
curl https://raw.githubusercontent.com/nk412/optparse/master/optparse.bash > optparse.bash

echo -e "${GREEN}  Installing up${NC}"
curl --create-dirs -o ~/.config/up/up.sh https://raw.githubusercontent.com/shannonmoeller/up/master/up.sh

has ansi > /dev/null
if [ $? -gt 0 ]; then
  echo -e "${GREEN}  Installing ansi${NC}"
  curl -OL git.io/ansi
  chmod 755 ansi
  sudo mv ansi /usr/local/bin/
fi

echo -e "${GREEN}  Installing loop${NC}"
snap install loop-rs --beta

npm_packages=""
npm_names=""

has how2 > /dev/null
if [ $? -gt 0 ]; then
  npm_packages="${npm_packages} how-2"
  npm_names="${npm_names} how2"
fi

has is > /dev/null
if [ $? -gt 0 ]; then
  npm_packages="${npm_packages} is.sh"
  npm_names="${npm_names} is"
fi

has mdlt > /dev/null
if [ $? -gt 0 ]; then
  npm_packages="${npm_packages} mdlt"
  npm_names="${npm_names} mdlt"
fi

has rename > /dev/null
if [ $? -gt 0 ]; then
  npm_packages="${npm_packages} rename-cli"
  npm_names="${npm_names} rename"
fi
echo -e "${GREEN}Installing the following packages using NPM: ${npm_names}${NC}"
sudo npm install -g ${npm_packages}

# Add in the project directories from github
# Atomsk
if [ -z $(which atomsk) ]; then
  cd ~/projects
  git clone git@github.com:pierrehirel/atomsk.git
  cd atomsk
  git remote add upstream git@github.com:pierrehirel/atomsk.git
  mr register
  cd src
  make atomsk
  sudo make install
else
  echo -e "${RED}Atomsk already installed at $(which atomsk)${NC}"
fi

# CXXOPTS
if [ ! -f /usr/local/include/cxxopts.hpp ]; then
  cd ~/projects
  git clone git@github.com:jarro2783/cxxopts.git
  cd cxxopts
  git remote add upstream git@github.com:jarro2783/cxxopts.git
  mr register
  sudo cp include/cxxopts.hpp /usr/local/include/
else
  echo -e "${RED}Cxxopts already installed${NC}"
fi

# has
if [ -z $(which has) ]; then
  cd ~/projects
  git clone git@github.com:kdabir/has.git
  cd has
  git remote add upstream git@github.com:kdabir/has.git
  mr register
  sudo make install
else
  echo -e "${RED}\"Has\" already installed at $(which has)${NC}"
fi

# LAMMPS
if [ ! -d ${HOME}/projects/lammps ]; then
  cd ~/projects
  git clone git@github.com:lammps/lammps.git
  cd lammps
  git remote add upstream git@github.com:lammps/lammps.git
  mr register
  echo -e "${GREEN}LAMMPS requires specific build directives that I cannot guess at - see the LAMMPS manual for more info.${NC}"
else
  echo -e "${RED}LAMMPS already installed - see LAMMPS manual for build info.${NC}"
fi

# MOOSE
if [ ! -d ${HOME}/projects/moose ]; then
  cd ~/projects
  git clone git@github.com:jarinfrench/moose.git
  cd moose
  git remote add upstream git@github.com:idaholab/moose.git
  mr register
  git fetch upstream
  git checkout devel
  git rebase upstream/devel
  ./scripts/update_and_rebuild_libmesh.
  cd modules/phase_field
  make -j24
  METHOD=dbg make -j24
else
  echo -e "${RED}MOOSE already installed${NC}"
fi

# Rtags
if [ -z $(which rdm) ]; then
  cd ~/projects
  git clone --recursive git@github.com:Andersbakken/rtags.git
  cd rtags
  git remote add upstream git@github.com:Andersbakken/rtags.git
  mr register
  mkdir build && cd build
  cmake -DCMAKE_EXPORT_COMPILE_COMMANDS=1 ..
  make
else
  echo -e "${RED}Rtags already installed - see the wiki for additional info.${NC}"
fi

unalias has
