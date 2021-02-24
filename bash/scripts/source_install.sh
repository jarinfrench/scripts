#! /bin/bash

shopt -s expand_aliases
box=$(uname -s)
RED='\033[0;31m'
GREEN='\e[0;32m'
NC='\033[0m'
alias has='curl -sL https://git.io/_has | HAS_ALLOW_UNSAFE=y bash -s'

if ! has curl snap git make wget; then
  echo -e "${RED}Some software will not be able to be installed unless the above is already installed${NC}"
fi

source_packages="optparse.bash up bd has ansi bashhub go"
if [[ "${box}" == "Linux" ]]; then
  for i in hr htop fzf hyperfine; do
    if ! has ${i} > /dev/null; then
      source_packages="${source_packages} ${i}"
    fi
  done
fi

echo -e "${GREEN}Installing the following packages from source: ${source_packages}${NC}"

for package in ${source_packages}; do
  echo -e "${GREEN}Installing ${package}${NC}"
  case "${package}" in
    optparse.bash)
      if ! has curl > /dev/null; then
        break
      else
        if [[ "${box}" == "Darwin" ]]; then
          echo "${GREEN}Note that optparse requires the GNU version of sed (install with brew install gnu-sed --with-default-names)${NC}"
        fi
        curl https://raw.githubusercontent.com/nk412/optparse/master/optparse.bash > optparse.bash
      fi
      ;;
    up)
      curl --create-dirs -o ~/.config/up/up.sh https://raw.githubusercontent.com/shannonmoeller/up/master/up.sh
      ;;
    hr)
      curl https://raw.githubusercontent.com/LuRsT/hr/master/hr | sudo tee -a hr > /dev/null
      sudo chmod +x hr
      sudo mv hr /usr/local/bin/
      ;;
    fzf)
      if [[ ! -d "${HOME}/.fzf" ]]; then
        git clone --depth 1 https://github.com/junegunn/fzf.git "${HOME}/.fzf"
        "${HOME}"/.fzf/install
      fi
      ;;
    bashhub)
      curl -OL https://bashhub.com/setup && bash setup
      rm setup
      ;;
    bd)
      if [[ -f "${HOME}/.bash_aliases" ]]; then
        # shellcheck source=/home/jarinf/.bash_aliases
        . "${HOME}/.bash_aliases"
      fi
      sudo wget --no-check-certificate -O /usr/local/bin/bd https://raw.github.com/vigneshwaranr/bd/master/bd
      sudo chmod +rx /usr/local/bin/bd
      # checks if the alias bd exists, and if not, adds the alias to the alias list.
      alias bd 2>/dev/null >/dev/null || (echo -e 'alias bd=". bd -si"' >> .bash_aliases && source .bash_aliases)
      ;;
    ansi)
      curl -OL git.io/ansi
      chmod 755 ansi
      sudo mv ansi /usr/local/bin
      ;;
    hyperfine)
      wget https://github.com/sharkdp/hyperfine/releases/download/v1.10.0/hyperfine-musl_1.10.0_amd64.deb
      sudo dpkg -i hyperfine-musl_1.10.0_amd64.deb
      rm hyperfine-musl_1.10.0_amd64.deb
      ;;
    has)
      curl -sL https://git.io/_has > ${HOME}/.local/bin/has
      ;;
    go)
      echo -e "${GREEN}Installing go1.14.6${NC}"
      if [[ "${box}" == "Linux" ]]; then
        go_package="https://golang.org/dl/go1.14.6.linux-amd64.tar.gz"
      elif [[ "${box}" == "Darwin" ]]; then
        go_package="https://golang.org/dl/go1.14.6.darwin-amd64.pkg"
      else
        echo -e "${RED}Unable to install go"
        continue
      fi
      wget -c "${go_package}" -O - | sudo tar -C /usr/local -xzf
      echo "Testing Go installation..."
      tmpfile=$(mktemp 2>/dev/null) || tmpfile=/tmp/go_test$$
      cat <<-EOF > "${tmpfile}"
        package main
        import "fmt"
        func main() {
          fmt.Printf("Hello world!\n")
        }
      EOF
      go build "${tmpfile}"
      tmp="${tmpfile%.go}"
      executable="${tmp##*/}"
      result="$(executable)"
      if [[ ! $(echo "${result}" | tr -d '[:space:]') == "Helloworld!" ]]; then
         echo "Go installation failed."
      fi
      rm "${tmpfile} ${executable}" 
        
    *)
      echo "Not implemented yet: ${package}"
      ;;
    esac
  done
