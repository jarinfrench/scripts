#! /bin/bash

shopt -s expand_aliases
RED='\033[0;31m'
GREEN='\e[0;32m'
NC='\033[0m'
snap_names=""

echo -e "${GREEN}Installing the following snaps: ${snap_names}${NC}"

for snp in ${snap_names}; do
  case ${snp} in
    *)
      echo -e "${RED}Not implemented yet: ${snp}${NC}"
      ;;
  esac
done
