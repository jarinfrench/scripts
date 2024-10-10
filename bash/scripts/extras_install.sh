#! /bin/bash

shopt -s expand_aliases
RED='\033[0;31m'
GREEN='\e[0;32m'
NC='\033[0m'
additional_scripts="anim2gif.sh gif.sh gif2anim.sh merge_gifs.sh"

echo "Installing the following scripts to ${HOME}/projects/scripts/bin:"
echo -e "${GREEN}${additional_scripts}${NC}"
