alias brls="git for-each-ref --count=30 --sort=-committerdate refs/heads/ --format='%(refname:short)'"
alias go="git checkout"
alias gs="git status"
alias gc="git commit"
alias ga="git add"
alias a="atom"
alias moose-update='cd ~/projects/moose/ && git fetch upstream && git checkout devel &&  git rebase upstream/devel; cd $OLDPWD'
alias start-moose='. ~/.moose_setup'
alias python-debug='python -m pdb'
alias python3-debug='python3 -m pdb'
alias fhere='find . -name'
alias bd=". bd -si" # See README for where to get this program

if [[ "$(uname -s)" = "Darwin" ]]; then
  # Mac specific aliases
  alias remove_ds_store='cd ~ && fhere "\.DS_Store" | xargs rm; cd ${OLDPWD}'
  alias alert='terminal-notifier -title "Terminal" -message "Done with task."; bash -c "printf '\a'"'
  alias moose-run="mpiexec -n 12 ~/projects/moose/modules/phase_field/phase_field-opt -i"
  alias moose-dbg='sudo lldb -- ~/projects/moose/modules/phase_field/phase_field-dbg -i'
  alias moose-make='cd ~/projects/moose/modules/phase-field && make -j24 && METHOD=dbg make -j24; cd ${OLDPWD}'
  alias lln="gls -lhtr  --time-style long-iso | gtac | gcat -n | gtac | gsed -s 's/^\s*\([0-9]*\)\s*\(.*\)/[\1]  \2 [\1]/'g && pwd"
elif [[ "$(uname -s)" = "Linux" ]]; then
  # Linux specific aliases
  alias moose-dbg='gdb -- ~/projects/moose/modules/phase_field/phase_field-dbg -i'
  alias moose-run="mpiexec -n 8 ~/projects/moose/modules/phase_field/phase_field-opt -i"
  alias moose-make='cd ~/projects/moose/modules/phase-field && make -j8 && METHOD=dbg make -j8; cd ${OLDPWD}'
  alias lln="ls -lhtr  --time-style long-iso | tac | cat -n | tac | sed -s 's/^\s*\([0-9]*\)\s*\(.*\)/[\1]  \2 [\1]/'g && pwd"
else
  echo "Unrecognized system: $(uname -s)"
  echo "No known $(uname -s)-specific aliases have been created."
fi
