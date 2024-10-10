#! /bin/bash
# shellcheck disable=SC2164

shopt -s expand_aliases
project_dir="${HOME}/projects"
RED='\033[0;31m'
REDBLD='\033[01;31m'
GREEN='\e[0;32m'
NC='\033[0m'

if [[ ! -d "${project_dir}" ]]; then
  echo -e "${RED}${project_dir} not found. Cannot continue installation${NC}"
  exit
fi
cur_dir=$(pwd)

git_repos="atomsk cxxopts has lammps moose rtags cod git_summary iprPy rng"
alias _has='curl -sL https://git.io/_has | HAS_ALLOW_UNSAFE=y bash -s'

echo "The following packages will be installed via git. The repos will be saved at ${project_dir}:"
echo -e "${GREEN}${git_repos}${NC}"

for repo in ${git_repos}; do
  case ${repo} in
    atomsk)
      if ! _has atomsk > /dev/null; then
        cd "${project_dir}"
        if [[ -d "${repo}" ]]; then
          echo -e "${RED}${repo} directory already exists${NC}"
        else
          echo -e "${GREEN}Cloning ${repo} github repository${NC}"
          git clone git@github.com:pierrehirel/atomsk.git
        fi
        cd atomsk
        git remote add upstream git@github.com:pierrehirel/atomsk.git
        if _has mr > /dev/null; then
          mr register
        fi
        cd src
        make atomsk && sudo make install
        cd "${project_dir}"
      else
        echo -e "${RED}atomsk already installed ($(which atomsk))${RED}"
      fi
      ;;
    cxxopts)
      if [[ ! -f /usr/local/include/cxxopts.hpp ]]; then
        cd "${project_dir}"
        if [[ -d "${repo}" ]]; then
          echo -e "${RED}${repo} directory already exists${NC}"
        else
          echo -e "${GREEN}Cloning ${repo} github repository${NC}"
          git clone git@github.com:jarro2783/cxxopts.git
        fi
        cd cxxopts
        git remote add upstream git@github.com:jarro2783/cxxopts.git
        if _has mr > /dev/null; then
          mr register
        fi
        echo -e "${GREEN}Installing header file to /usr/local/include/cxxopts.hpp${NC}"
        sudo cp include/cxxopts.hpp /usr/local/include/
        cd "${project_dir}"
      else
        echo -e "${RED}cxxopts already installed at /usr/local/include/cxxopts.hpp${NC}"
      fi
      ;;
    has)
      if ! _has has > /dev/null; then
        cd "${project_dir}"
        if [[ -d "${repo}" ]]; then
          echo -e "${RED}${repo} directory already exists${NC}"
        else
          echo -e "${GREEN}Cloning ${repo} github repository${NC}"
          git clone git@github.com:kdabir/has.git
        fi
        cd has
        git remote add upstream git@github.com:kdabir/has.git
        if _has mr > /dev/null; then
          mr register
        fi
        sudo make install
        cd "${project_dir}"
      else
        echo -e "${RED}\"has\" already installed ($(which has))${NC}"
      fi
      ;;
    lammps)
      if [[ -d "${project_dir}/lammps" ]]; then
        cd "${project_dir}"
        git clone git@github.com:lammps/lammps.git

        cd lammps
        git remote add upstream git@github.com:lammps/lammps.git
        if _has mr > /dev/null; then
          mr register
        fi
        echo -e "${GREEN}See LAMMPS manual for specific build information${NC}"
        cd "${project_dir}"
      else
        echo -e "${RED}LAMMPS already installed - see LAMMPS manual for build information${NC}"
      fi
      ;;
    moose)
      if [[ ! -d "${project_dir}/moose" ]]; then
        cd "${project_dir}"
        echo -e "${GREEN}Cloning ${repo} github repository${NC}"
        git clone git@github.com:jarinfrench/moose.git

        cd moose
        git remote add upstream git@github.com:idaholab/moose.git
        if _has mr > /dev/null; then
          mr register
        fi
        git fetch upstream
        git checkout devel
        git rebase upstream/devel
        echo -e "${GREEN}Building libmesh${NC}"
        ./scripts/update_and_rebuild_libmesh.sh
        cd modules/phase_field
        echo -e "${GREEN}Building phase_field module${NC}"
        make -j "$(grep -c processor /proc/cpuinfo)"
        METHOD=dbg make -j "$(grep -c processor /proc/cpuinfo)"
        cd "${project_dir}"
      else
        echo -e "${RED}MOOSE already installed${NC}"
      fi
      ;;
    rtags)
      if ! _has rdm > /dev/null; then
        cd "${project_dir}"
        if [[ -d "${repo}" ]]; then
          echo -e "${RED}${repo} directory already exists${NC}"
        else
          echo -e "${GREEN}Cloning ${repo} github repository${NC}"
          git clone --recursive git@github.com:Andersbakken/rtags.git
        fi
        cd rtags
        git remote add upstream git@github.com:Andersbakken/rtags.git
        if _has mr > /dev/null; then
          mr register
        fi
        mkdir build; cd build
        echo -e "${GREEN}Building rtags${NC}"
        cmake -DCMAKE_EXPORT_COMPILE_COMMANDS=1 ..
        make
        cd "${project_dir}"
      else
        echo -e "${RED}rtags already installed - see the wiki for additional info${NC}"
      fi
      ;;
    cod)
      if ! _has cod > /dev/null; then
        cd "${project_dir}"
        if [[ -d "${repo}" ]]; then
          echo -e "${RED}${repo} directory already exists${NC}"
        else
          echo -e "${GREEN}Cloning ${repo} github repository${NC}"
          git clone https://github.com/dim-an/cod.git
        fi
        cd cod
        if _has go > /dev/null; then
          go build
          echo -e "${GREEN}Installing cod to /usr/local/bin${NC}"
          sudo cp cod /usr/local/bin
        else
          echo -e "${RED}Go required to build cod${NC}"
        fi
        if _has mr > /dev/null; then
          mr register
        fi
        git remote add upstream git@github.com:dim-an/cod.git
        if ! grep 'source <(cod init $$ bash)' "${HOME}/.bashrc"; then
          echo -e "${GREEN}" 'Adding the line "source <(cod init $$ bash)"' "to ${HOME}/.bashrc${NC}"
          echo 'source <(cod init $$ bash)' >> "${HOME}/.bashrc"
        fi
        cd "${project_dir}"
      else
        echo -e "${RED}cod already installed${NC}"
      fi
      ;;
    git_summary)
      if ! _has git-summary; then
        cd "${project_dir}"
        if [[ -d "git-summary" ]]; then
          echo -e "${RED}git-summary directory already exists${NC}"
        else
          echo -e "${GREEN}Cloning ${repo} github repository${NC}"
          git clone git@github.com:MirkoLedda/git-summary.git
        fi
        cd git-summary
        if _has mr > /dev/null; then
          mr register
        fi
        git remote add upstream git@github.com:MirkoLedda/git-summary.git
        echo -e "${GREEN}Installing git-summary to /usr/local/bin${NC}"
        sudo cp git-summary /usr/local/bin
        cd "${project_dir}"
      else
        echo -e "${RED}git-summary already installed${NC}"
      fi
      ;;
    iprPy)
      if [[ ! -d "${project_dir}/iprPy" ]]; then
        cd "${project_dir}"
        git clone git@github.com:jarinfrench/iprPy.git

        cd iprPy
        git remote add upstream git@github.com:lmhale99/iprPy.git
        if _has mr > /dev/null; then
          mr register
        fi
        # python setup.py develop --user
        # cd bin
        echo -e "${REDBLD}NOTE that as of 19 June 2020, iprPy still uses Python 2, and thus will not run!${NC}"
        cd "${project_dir}"
      else
        echo -e "${RED}ipyPr already downloaded${NC}"
      fi
      ;;
    rng)
      if [[ ! -d "${project_dir}/rng" ]]; then
        cd "${project_dir}"
        git clone https://github.com/nickolasburr/rng.git
        cd rng
        make && sudo make install
        cd "${project_dir}"
      else
        echo -e "${RED}rng already installed${NC}"
      fi
      ;;
    *)
      echo -e "${RED}${repo} installation instructions not written${NC}"
      ;;
  esac
done

cd "${cur_dir}"
unalias _has
