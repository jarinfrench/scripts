#! /bin/bash

box=$(uname -s)

if ! command -v curl &> /dev/null; then
  echo "curl is required for this install, but is not found. Installing curl..."
  if command -v brew &> /dev/null; then
    brew install curl
  else
    echo "Homebrew not found on ${box} system, unable to install curl. Exiting..."
    exit 1
  fi
fi

shopt -s expand_aliases
RED='\033[0;31m'
GREEN='\e[0;32m'
NC='\033[0m'
alias has='curl -sL https://git.io/_has | HAS_ALLOW_UNSAFE=y bash -s'

if [[ "${box}" != "Darwin" ]]; then
  exit
fi

if ! has brew; then
  echo "${RED}Homebrew required for installation${NC}"
  exit
fi
echo -e "Mac system detected, installing packages using Homebrew"

brew_packages=""
for i in fd hr mr jrnl task taskd tasksh rng cloc fzf hyperfine coreutils; do
  if ! has ${i} > /dev/null; then
    if [ "${i}" == "rng" ]; then
      tap_brew=y
    fi
    brew_packages="${brew_packages} ${i}"
  fi
done

echo -e "The following packages will be installed from Homebrew: ${GREEN}${brew_packages}${NC}"
if [[ "${tap_brew}" ]]; then
  brew tap nickolasburr/pfa
fi
brew install ${brew_packages}
