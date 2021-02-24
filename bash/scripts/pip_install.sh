#! /bin/bash

box=$(uname -s)

shopt -s expand_aliases
alias has='curl -sL https://git.io/_has | HAS_ALLOW_UNSAFE=y bash -s'

if ! has pip3 > /dev/null; then
  echo "pip3 required for installing pip packages"
  exit 1
fi

pip_packages=""
if [[ "${box}" == "Linux" ]]; then
  pip_packages="${pip_packages} jrnl"
fi

################################################################################
# jrnl - command line journal keeping software. See https://pypi.org/project/jrnl/
# argparse - an argument parsing library. See https://docs.python.org/3/library/argparse.html
# tqdm - Progress bar for python and CLI. See https://github.com/tqdm/tqdm
# colorama - ANSI character escape sequences made easy. See https://pypi.org/project/colorama/
# natsort - a natural sorting library. See https://pypi.org/project/natsort/
# prospector - Analyzes python code. See https://github.com/PyCQA/prospector
# vulture - Finds dead (unused) python code. See https://github.com/jendrikseipp/vulture
# pylint - Python linter. See https://www.pylint.org/
# pytest - Python testing framework. See https://docs.pytest.org/en/latest/
# coverage - Examine coverage of tests. See https://pypi.org/project/coverage/
# MDAnalysis - Library to examine MD trajectories. See https://www.mdanalysis.org/
# atomman - Atomistic manipulation toolkit. See https://github.com/usnistgov/atomman/
################################################################################
for i in argparse tqdm colorama natsort prospector vulture pylint pytest coverage MDAnalysis atomman; do
  pip_packages="${pip_packages} ${i}"
done

echo "Installing pip packages: ${pip_packages}"
pip install --user -U "${pip_packages}" # also upgrades any existing packages
