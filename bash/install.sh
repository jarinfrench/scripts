#! /bin/bash

function ctrl_c() {
  echo -e "Cleaning up..."
  if [[ -f ${HOME}/.bash_aliases.bak ]]; then
    mv ~/.bash_aliases.bak ~/.bash_aliases
  fi
  if [[ -f ${HOME}/.bash_functions.bak ]]; then
    mv ~/.bash_functions.bak ~/.bash_functions
  fi
  if [[ -f ${HOME}/.bash_variables.bak ]]; then
    mv ~/.bash_variables.bak ~/.bash_variables
  fi
  if [[ -f ${HOME}/.bashrc.bak ]]; then
    mv ~/.bashrc.bak ~/.bashrc
  fi
  if [[ -f ${HOME}/.bash_profile.bak ]]; then
    mv ~/.bash_profile.bak ~/.bash_profile
  fi
  if [[ -f ${HOME}/.pdbrc.bak ]]; then
    mv ~/.pdbrc.bak ~/.pdbrc
  fi
  if [[ -f ${HOME}/.jrnl_config.bak ]]; then
    mv ~/.jrnl_config.bak ~/.jrnl_config
  fi
  if [[ -f ${HOME}/.mrconfig.bak ]]; then
    mv ~/.mrconfig.bak ~/.mrconfig
  fi
  if [[ -f ${HOME}/.taskrc.bak ]]; then
    mv ~/.taskrc.bak ~/.taskrc
  fi

  exit 1
}

GREEN='\033[0;32m'
RED='\033[0;31m'
NC='\033[0m'

trap ctrl_c INT

INSTALL_DIR="${HOME}/projects/scripts/bash"
cd ${INSTALL_DIR}

# -b uses a simple backup for the destination file (if it exists), and -S changes the suffix from ~ to .bak (in this case)
for i in bash_aliases bash_functions bash_variables pdbrc jrnl_config mrconfig taskrc vimrc; do
  if [ -f .${i}.bak ]; then
    echo ".${i} already has a backup file (.${i}.bak) - unable to link file."
  else
    ln -sbS .bak ${INSTALL_DIR}/.${i} ~/.${i}
  fi
done

ln -sf ${INSTALL_DIR}/.alias_completion.sh ~/.alias_completion.sh

# Store symlinks to python programs in bin
for i in $(find ~/projects/scripts/python/ -type f -executable); do
  file=$(basename ${i})
  ln -sf ${i} ~/projects/scripts/bin/${file}
done

UNAME=$(uname -s)
if [[ "${UNAME}" == "Darwin" ]]; then
  echo -e "${GREEN}If you would like your current .bash_profile overwritten, execute the following command:"
  echo -e "ln -sbS .bak ${INSTALL_DIR}/.bash_profile ~/.bash_profile"
  echo -e "${NC}"

  echo -e "${GREEN}If you would like your current .bashrc overwritten, execute the following command:"
  echo -e "ln -sbS .bak ${INSTALL_DIR}/.bashrc_mac ~/.bashrc"
  echo -e "${NC}"

  # required programs:
  curl -sL https://git.io/_has | bash -s snap curl brew wget npm git
  if [ $? -gt 0 ]; then
    echo -e "${RED}Make sure you have the above software installed!${NC}"
    exit 2
  fi
elif [[ "${UNAME}" == "Linux" ]]; then
  echo -e "${GREEN}If you would like your current .bashrc overwritten, execute the following command:"
  echo -e "ln -sbS .bak ${INSTALL_DIR}/.bashrc_linux ~/.bashrc"
  echo -e "${NC}"

  curl -sL https://git.io/_has | bash -s snap curl wget npm git pip3
  if [ $? -gt 0 ]; then
    echo -e "${RED}Make sure you have the above software installed!${NC}"
    exit 2
  fi
fi

if [[ "${UNAME}" == "Darwin" ]]; then
  echo -e "Mac system detected, installing packages using Homebrew."
  # install using Darwin methods
  # Note that these methods have been tested on Sierra, and found to work!
  echo -e "Installing the following packages from Homebrew: fd hr mr jrnl task taskd tasksh rng cloc"
  brew tap nickolasburr/pfa
  brew install fd hr mr jrnl task taskd tasksh rng cloc htop

  echo -e "Installing the following packages from source: bd has optparse.bash up ansi"

elif [[ "${UNAME}" == "Linux" ]]; then
  echo -e "${GREEN}Linux system detected, installing packages via apt."
  # install using UBUNTU methods (may need to change this later, but it should work for now)
  echo -e "Installing the following packages from repositories: myrepos taskwarrior cloc htop parallel${NC}"
  sudo apt install myrepos taskwarrior cloc htop parallel

  echo -e "  ${GREEN}Configuring taskwarrior... please type yes when prompted${NC}"
  ln -sf ${INSTALL_DIR}/.task ${HOME}/.task
  task config taskd.server freecinc.com:53589
  task config taskd.key ~/.task/freecinc_b77224d4.key.pem
  task config taskd.certificate ~/.task/freecinc_b77224d4.cert.pem
  task config taskd.ca ~/.task/freecinc_b77224d4.ca.pem
  task config taskd.credentials -- 'FreeCinc/freecinc_b77224d4/cfdcc1e7-359c-4291-b7cc-34f91ad4d8b6'

  echo -e "${GREEN}Installing the following package using pip: jrnl${NC}"
  pip3 install --user jrnl

  echo -e "${GREEN}Installing the following packages from source: bd fd has hr up htop ansi optparse.bash"
  echo -e "\tInstalling fd${NC}"
  curl https://github.com/sharkdp/fd/releases/download/v7.3.0/fd-musl_7.3.0_amd64.deb > fd-musl_7.3.0_amd64.deb
  sudo dpkg -i fd-musl_7.3.0_amd64.deb
  rm fd-musl_7.3.0_amd64.deb

  echo -e "${GREEN}\tInstalling hr${NC}"
  sudo curl https://raw.githubusercontent.com/LuRsT/hr/master/hr > hr
  sudo chmod +x hr
  sudo mv hr /usr/local/bin/

  if [ ! -d ${HOME}/projects/scripts/bash/rng ]; then
    echo -e "${GREEN}\tInstalling rng${NC}"
    git clone https://github.com/nickolasburr/rng.git
    cd rng
    make
    sudo make install
    cd ../
  fi
fi

. .bash_aliases # required to check if bd is aliased

echo -e "${GREEN}\tInstalling bd${NC}"
sudo wget --no-check-certificate -O /usr/local/bin/bd https://raw.github.com/vigneshwaranr/bd/master/bd
sudo chmod +rx /usr/local/bin/bd
alias bd 2>/dev/null >/dev/null || (echo -e 'alias bd=". bd -si"' >> .bash_aliases && source .bash_aliases) # checks if the alias bd exists, and if not, adds the alias to the alias list.

echo -e "${GREEN}\tInstalling optparse.bash"
echo -e "\tNote that optparse requires the GNU version of sed (for Mac - install by brew install gnu-sed (--with-default-names, if you don't want to alias sed))${NC}"
curl https://raw.githubusercontent.com/nk412/optparse/master/optparse.bash > optparse.bash

echo -e "${GREEN}\tInstalling up${NC}"
curl --create-dirs -o ~/.config/up/up.sh https://raw.githubusercontent.com/shannonmoeller/up/master/up.sh

echo -e "${GREEN}\tInstalling ansi${NC}"
curl -OL git.io/ansi
chmod 755 ansi
sudo mv ansi /usr/local/bin/

echo -e "${GREEN}Installing loop${NC}"
snap install loop-rs --beta

echo -e "${GREEN}Installing the following packages using NPM: how2 is${NC}"
sudo npm install -g how-2 is.sh

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
