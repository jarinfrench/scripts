#! /bin/bash

INSTALL_DIR="~/projects/scripts/bash"
cd "${INSTALL_DIR}"

ln -s .bash_aliases ~/.bash_aliases
ln -s .bash_functions ~/.bash_functions
ln -s .bash_variables ~/.bash_variables
ln -s .pdbrc ~/.pdbrc
ln -s .jrnl_config ~/.jrnl_config
ln -s .mrconfig ~/.mrconfig
ln -s .taskrc ~/.taskrc

# Store symlinks to python programs in bin
# Note that there is probably a better way to do this involving searching for the shebang in the python directory
ln -s ~/projects/scripts/python/average_data.py ~/projects/scripts/bin/average_data.py
ln -s ~/projects/scripts/python/calculate_force_and_velocity.py ~/projects/scripts/bin/calculate_force_and_velocity.py
ln -s ~/projects/scripts/python/calculate_lattice_param.py ~/projects/scripts/bin/calculate_lattice_param.py
ln -s ~/projects/scripts/python/calculate_statistics.py ~/projects/scripts/bin/calculate_statistics.py
ln -s ~/projects/scripts/python/convertOmatToEuler.py ~/projects/scripts/bin/convertOmatToEuler.py
ln -s ~/projects/scripts/python/Euler_angle_generator.py ~/projects/scripts/bin/Euler_angle_generator.py
ln -s ~/projects/scripts/python/generate_polycrystal_data.py ~/projects/scripts/bin/generate_polycrystal_data.py
ln -s ~/projects/scripts/python/orientation_matrix.py ~/projects/scripts/bin/orientation_matrix.py
ln -s ~/projects/scripts/python/parse_grain_growth_output.py ~/projects/scripts/bin/parse_grain_growth_output.py
ln -s ~/projects/scripts/python/plot_data.py ~/projects/scripts/bin/plot_data.py
ln -s ~/projects/scripts/python/plot_grain_sizes.py ~/projects/scripts/bin/plot_grain_sizes.py
ln -s ~/projects/scripts/python/plot_LAMMPS_data.py ~/projects/scripts/bin/plot_LAMMPS_data.py
ln -s ~/projects/scripts/python/rotation_matrix.py ~/projects/scripts/bin/rotation_matrix.py
ln -s ~/projects/scripts/python/uo2_lattice_param.py ~/projects/scripts/bin/uo2_lattice_param.py

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
  echo "Mac system detected, installing packages using Homebrew."
  # install using Darwin methods
  # Note that these methods have been tested on Sierra, and found to work!
  echo "Installing the following packages from Homebrew: fd hr mr jrnl task taskd tasksh rng cloc"
  brew tap nickolasburr/pfa
  brew install fd hr mr jrnl task taskd tasksh rng cloc htop

  echo "Installing the following packages from source: bd has optparse.bash up ansi"
  echo -e "\tInstalling bd"
  wget --no-check-certificate -O /usr/local/bin/bd https://raw.github.com/vigneshwaranr/bd/master/bd
  chmod +rx /usr/local/bin/bd
  alias bd 2>/dev/null >/dev/null || (echo 'alias bd=". bd -si"' >> .bash_aliases && source .bash_aliases) # checks if the alias bd exists, and if not, adds the alias to the alias list.

elif [[ "${UNAME}" == "Linux" ]]; then
  echo "Linux system detected, installing packages via apt."
  # install using UBUNTU methods (may need to change this later, but it should work for now)
  echo "Installing the following packages from repositories: bd myrepos taskwarrior cloc htop"
  sudo apt install bd myrepos taskwarrior cloc htop

  echo "Installing the following packages from source: fd has hr up htop ansi optparse.bash"
  echo -e "\tInstalling fd"
  wget https://github.com/sharkdp/fd/releases/download/v7.3.0/fd-musl_7.3.0_amd64.deb
  sudo dpkg -i fd-musl_7.3.0_amd64.deb

  echo -e "\tInstalling hr"
  curl https://raw.githubusercontent.com/LuRsT/hr/master/hr > ~/bin/hr
  (Examine ~/bin/hr)
  chmod +x ~/bin/hr

  echo -e "\tInstalling rng"
  git clone https://github.com/nickolasburr/rng.git
  cd rng
  make
  sudo make install
  cd ../

  echo "Installing the following package using pip: jrnl "
  pip3 install --user jrnl
fi

echo -e "\tInstalling optparse.bash"
echo -e "\tNote that optparse requires the GNU version of sed (for Mac - install by brew install gnu-sed (--with-default-names, if you don't want to alias sed))"
wget https://raw.githubusercontent.com/nk412/optparse/master/optparse.bash

echo -e "\tInstalling up"
curl --create-dirs -o ~/.config/up/up.sh https://raw.githubusercontent.com/shannonmoeller/up/master/up.sh

echo -e "\tInstalling ansi"
curl -OL git.io/ansi
chmod 755 ansi
sudo mv ansi /usr/local/bin/

echo "Installing the following packages using NPM: how2 is"
npm install -g how-2 is.sh

# Add in the project directories from github
# Atomsk
cd ~/projects
git clone git@github.com:pierrehirel/atomsk.git
cd atomsk
git remote add upstream git@github.com:pierrehirel/atomsk.git
mr register
cd src
make atomsk
sudo make install

# CXXOPTS
cd ~/projects
git clone git@github.com:jarro2783/cxxopts.git
cd cxxopts
git remote add upstream git@github.com:jarro2783/cxxopts.git
mr register
sudo cp include/cxxopts.hpp /usr/local/include/

# has
cd ~/projects
git clone git@github.com:kdabir/has.git
cd has
git remote add upstream git@github.com:kdabir/has.git
mr register
sudo make install

# LAMMPS
cd ~/projects
git clone git@github.com:lammps/lammps.git
cd lammps
git remote add upstream git@github.com:lammps/lammps.git
mr register
echo "LAMMPS requires specific build directives that I cannot guess at - see the LAMMPS manual for more info."

# MOOSE
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

# Rtags
cd ~/projects
git clone --recursive git@github.com:Andersbakken/rtags.git
cd rtags
git remote add upstream git@github.com:Andersbakken/rtags.git
mr register
mkdir build && cd build
cmake -DCMAKE_EXPORT_COMPILE_COMMANDS=1 ..
make

# Uncrustify (for atom-beautify c++)
cd ~/projects
git clone git@github.com:uncrustify/uncrustify.git
cd uncrustify
git remote add upstream git@github.com:uncrustify/uncrustify.git
mr register
mkdir build && cd build
cmake ..
make
sudo make install
