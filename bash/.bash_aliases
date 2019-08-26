alias a="atom"
alias bd=". bd -si" # See README for where to get this program
alias brls="git for-each-ref --count=30 --sort=-committerdate refs/heads/ --format='%(refname:short)'"
alias chmod='chmod u+x' # gives the file owner (u) the additional (+) execute bit (x), making the specified file executable
alias fhere='find . -name'
alias gbstudio="appletviewer https://staff.aist.go.jp/h.ogawa/GBstudio/gbs3/GBs_PJ.html"
alias ga="git add"
alias gc="git commit"
alias go="git checkout"
alias gs="git status"
alias lh='ls -lisAd .[^.]*' # list just the hidden files in a directory
alias moose-commit='cd ~/projects/moose/ && go devel && git push origin devel && cd ${OLDPWD}'
alias moose-update='cd ~/projects/moose/ && git fetch upstream && git checkout devel &&  git rebase upstream/devel; cd $OLDPWD'
alias plot='plot_data.py'
alias python-debug='python -m pdb'
alias python3-debug='python3 -m pdb'
alias start-moose='. ~/.moose_setup'
alias 7za='7z a -m0=lzma -mx=9' # arguments needed are archive name, and files to compress

if [[ "$(uname -s)" = "Darwin" ]]; then
  # Mac specific aliases
  alias alert='terminal-notifier -title "Terminal" -message "Done with task."; bash -c "printf '\a'"'
  alias ctags="`brew --prefix`/bin/ctags"
  alias lln="gls -lhtr  --time-style long-iso | gtac | gcat -n | gtac | gsed -s 's/^\s*\([0-9]*\)\s*\(.*\)/[\1]  \2 [\1]/'g && pwd"
  alias moose-debug='sudo lldb -- ~/projects/moose/modules/phase_field/phase_field-dbg -i'
  alias moose-make='cd ~/projects/moose/modules/phase_field && make -j24 && METHOD=dbg make -j24; cd ${OLDPWD}'
  alias moose-run="mpiexec -n 12 ~/projects/moose/modules/phase_field/phase_field-opt -i"
  alias remove_ds_store='cd ~ && fhere "\.DS_Store" | xargs rm; cd ${OLDPWD}'
  alias sed="gsed"
elif [[ "$(uname -s)" = "Linux" ]]; then
  # Linux specific aliases
  alias lln="ls -lhtr  --time-style long-iso | tac | cat -n | tac | sed -s 's/^\s*\([0-9]*\)\s*\(.*\)/[\1]  \2 [\1]/'g && pwd"
  alias moose-debug='gdb -- ~/projects/moose/modules/phase_field/phase_field-dbg -i'
  alias moose-make='cd ~/projects/moose/modules/phase_field && make -j8 && METHOD=dbg make -j8; cd ${OLDPWD}'
  alias moose-run="mpiexec -n 8 ~/projects/moose/modules/phase_field/phase_field-opt -i"
  alias open='xdg-open'
  alias rm_junk='find . \( -name "*.[eo][1-9]*" -or -name "slurm-[1-9]*.out" \) | xargs rm'
  alias update='sudo apt-get update && sudo apt-get upgrade && sudo apt-get autoremove'
else
  echo "Unrecognized system: $(uname -s)"
  echo "No known $(uname -s)-specific aliases have been created."
fi
