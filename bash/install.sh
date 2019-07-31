#! /bin/bash

cd ~/projects/scripts/bash

ln -s .bash_aliases ~/.bash_aliases
ln -s .bash_functions ~/.bash_functions
ln -s .bash_variables ~/.bash_variables
ln -s .pdbrc ~/.pdbrc
ln -s .jrnl_config ~/.jrnl_config
ln -s .mrconfig ~/.mrconfig
ln -s .taskrc ~/.taskrc

local UNAME=$(uname -s)
if [[ "${UNAME}" == "Darwin" ]]; then
  echo "If you would like your current .bash_profile overwritten, execute the following command:"
  echo "ln -s .bash_profile ~/.bash_profile"
  echo ""
fi

echo "If you would like your current .bashrc overwritten, execute the following command:"
echo "ln -s .bashrc ~/.bashrc"
echo ""

if [[ "${UNAME}" == "Darwin" ]]; then
  # install using Darwin methods
  # Note that these methods have been tested on Sierra, and found to work!
  echo "Installing the following packages from Homebrew: fd hr mr jrnl task taskd tasksh rng cloc"
  brew tap nickolasburr/pfa
  brew install fd hr mr jrnl task taskd tasksh rng cloc htop

  echo "Installing the following packages from NPM: how2 is"
  npm install -g how-2 is.sh

  echo "Installing the following packages from source: bd has optparse.bash up ansi"
  echo -e "\tInstalling bd"
  wget --no-check-certificate -O /usr/local/bin/bd https://raw.github.com/vigneshwaranr/bd/master/bd
  chmod +rx /usr/local/bin/bd
  alias bd 2>/dev/null >/dev/null || (echo 'alias bd=". bd -si"' >> .bash_aliases && source .bash_aliases) # checks if the alias bd exists, and if not, adds the alias to the alias list.

  echo -e "\tInstalling has"
  cd ~/projects && git clone https://github.com/kdabir/has.git && cd has && sudo make install; cd ~/projects/scripts/bash

  echo -e "\tInstalling optparse.bash"
  echo -e "\tNote that optparse requires the GNU version of sed (install by brew install gnu-sed (--with-default-names, if you don't want to alias sed))"
  wget https://raw.githubusercontent.com/nk412/optparse/master/optparse.bash

  echo -e "\tInstalling up"
  curl --create-dirs -o ~/.config/up/up.sh https://raw.githubusercontent.com/shannonmoeller/up/master/up.sh

  echo -e "\tInstalling ansi"
  curl -OL git.io/ansi
  chmod 755 ansi
  sudo mv ansi /usr/local/bin/
elif [[ "${UNAME}" == "Linux" ]]; then
  # install using UBUNTU methods (may need to change this later, but it should work for now)
  echo "Installing the following packages from repositories: bd myrepos taskwarrior cloc htop"
  sudo apt-get install bd myrepos taskwarrior cloc htop

  echo "Installing the following packages from source: fd has hr up htop ansi optparse.bash"
  echo -e "\tInstalling fd"
  wget https://github.com/sharkdp/fd/releases/download/v7.3.0/fd-musl_7.3.0_amd64.deb
  sudo dpkg -i fd-musl_7.3.0_amd64.deb

  echo -e "\tInstalling has"
  cd ~/projects && git clone https://github.com/kdabir/has.git && cd has && sudo make install; cd ~/projects/scripts/bash

  echo -e "\tInstalling hr"
  curl https://raw.githubusercontent.com/LuRsT/hr/master/hr > ~/bin/hr
  (Examine ~/bin/hr)
  chmod +x ~/bin/hr

  echo -e "\tInstalling up"
  curl --create-dirs -o ~/.config/up/up.sh https://raw.githubusercontent.com/shannonmoeller/up/master/up.sh

  echo -e "\tInstalling ansi"
  curl -OL git.io/ansi
  chmod 755 ansi
  sudo mv ansi /usr/local/bin/

  echo -e "\tInstalling optparse.bash"
  wget https://raw.githubusercontent.com/nk412/optparse/master/optparse.bash

  echo -e "\tInstalling rng"
  git clone https://github.com/nickolasburr/rng.git
  cd rng
  make
  sudo make install
  cd ../

  echo "Installing the following packages using NPM: how2 is"
  npm install -g how-2
  npm install -g is.sh

  echo "Installing the following package using pip: jrnl "
  pip3 install --user jrnl

fi
