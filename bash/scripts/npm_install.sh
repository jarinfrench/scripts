#! /bin/bash

shopt -s expand_aliases
alias has='curl -sL https://git.io/_has | HAS_ALLOW_UNSAFE=y bash -s'
npm_packages=""
npm_names=""

if ! has how2 > /dev/null; then
  npm_packages="${npm_packages} how-2"
  npm_names="${npm_names} how2"
fi

if ! has is > /dev/null; then
  npm_packages="${npm_packages} is.sh"
  npm_names="${npm_names} is"
fi

if ! has mdlt > /dev/null; then
  npm_packages="${npm_packages} mdlt"
  npm_names="${npm_names} mdlt"
fi

if ! has rename > /dev/null; then
  npm_packages="${npm_packages} rename-cli"
  npm_names="${npm_names} rename"
fi

if ! has fd > /dev/null; then
  npm_packages="${npm_packages} fd-find"
  npm_names="${npm_names} fd"
fi

if [[ ! -z "${npm_packages}" ]]; then
  echo "Installing the following packages using npm: ${npm_names}"
  sudo npm install -g ${npm_packages} --unsafe-perm
fi
