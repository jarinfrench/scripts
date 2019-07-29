#! /bin/bash

ln -s .bash_aliases ~/.bash_aliases
ln -s .bash_functions ~/.bash_functions
ln -s .bash_variables ~/.bash_variables
ln -s .pdbrc ~/.pdbrc

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
  # install using Darwin methods... homebrew?

  # wget --no-check-certificate -O /usr/local/bin/bd https://raw.github.com/vigneshwaranr/bd/master/bd
  # chmod +rx /usr/local/bin/bd
  # echo 'alias bd=". bd -si"' >> ~/.bashrc
  # source ~/.bashrc
  # brew install fd
  # brew install hr
  # brew install mr
  # brew install jrnl
  # brew install task taskd tasksh
  # brew tap nickolasburr/pfa
  # brew install rng
  # brew install cloc
  # npm install -g how-2
  # npm install -g is.sh
  # wget https://raw.githubusercontent.com/nk412/optparse/master/optparse.bash
  # curl --create-dirs -o ~/.config/up/up.sh https://raw.githubusercontent.com/shannonmoeller/up/master/up.sh
  #
  # wget https://hisham.hm/htop/releases/2.2.0/htop-2.2.0.tar.gz
  # tar -xzvf htop-2.2.0.tar.gz
  # mv htop-2.2.0.tar.gz htop-2.2.0
  # cd htop-2.2.0
  # ./configure && make && sudo make install
  # cd ../
  # mv htop-2.2.0 ~/
  #
  # curl -OL git.io/ansi
  # chmod 755 ansi
  # sudo mv ansi /usr/local/bin/
elif [[ "${UNAME}" == "Linux" ]]; then
  # install using UBUNTU methods (may need to change this later, but it should work for now)
  echo "Attempting to install the following packages from repositories: bd myrepos taskwarrior cloc"
  sudo apt-get install bd myrepos taskwarrior cloc

  echo "Installing the following packages from source: fd hr up htop ansi optparse.bash"
  echo -e "\tInstalling fd"
  wget https://github.com/sharkdp/fd/releases/download/v7.3.0/fd-musl_7.3.0_amd64.deb
  sudo dpkg -i fd-musl_7.3.0_amd64.deb

  echo -e "\tInstalling hr"
  curl https://raw.githubusercontent.com/LuRsT/hr/master/hr > ~/bin/hr
  (Examine ~/bin/hr)
  chmod +x ~/bin/hr

  echo -e "\tInstalling up"
  curl --create-dirs -o ~/.config/up/up.sh https://raw.githubusercontent.com/shannonmoeller/up/master/up.sh

  echo -e "\tInstalling htop"
  wget https://hisham.hm/htop/releases/2.2.0/htop-2.2.0.tar.gz
  tar -xzvf htop-2.2.0.tar.gz
  mv htop-2.2.0.tar.gz htop-2.2.0
  cd htop-2.2.0
  ./configure && make && sudo make install
  cd ../
  mv htop-2.2.0 ~/

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
