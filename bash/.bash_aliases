alias brls="git for-each-ref --count=30 --sort=-committerdate refs/heads/ --format='%(refname:short)'"
alias go="git checkout"
alias gs="git status"
alias gc="git commit"
alias ga="git add"
alias a="atom"
alias moose-run="mpiexec -n 12 ~/projects/moose/modules/phase_field/phase_field-opt -i"
alias moose-dbg='sudo lldb -- ~/projects/moose/modules/phase_field/phase_field-dbg -i'
alias moose-update='cd ~/projects/moose/ && git fetch upstream && git checkout devel &&  git rebase upstream/devel && cd $OLDPWD'
alias moose-make='cd ~/projects/moose/modules/phase-field && make -j24 && METHOD=dbg make -j24 && cd ${OLDPWD}'
alias start-moose='. ~/.moose_setup'
alias alert='terminal-notifier -title "Terminal" -message "Done with task."; bash -c "printf '\a'"'
alias python-debug='python -m pdb'
alias python3-debug='python3 -m pdb'
alias fhere='find . -name'
alias lln="gls -lhtr  --time-style long-iso | gtac | gcat -n | gtac | gsed -s 's/^\s*\([0-9]*\)\s*\(.*\)/[\1]  \2 [\1]/'g && pwd"
alias bd=". bd -si" # See README for where to get this program
