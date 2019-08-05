# make bash append rather than overwrite the history on disk:
shopt -s histappend
export HISTCONTROL=ignoredups:erasedups
export HISTSIZE=100000
export HISTFILESIZE=100000

# Whenever displaying the prompt, write the previous line to disk
export PROMPT_COMMAND="${PROMPT_COMMAND:+$PROMPT_COMMAND ;} history -a"

# When changing directories, small types can be ignored by bash
shopt -s cdspell

# source any user-defined bash aliases in .bash_aliases
if [ -f ~/.bash_aliases ]; then
  . ~/.bash_aliases
fi

source ~/projects/scripts/bash/git_completion.sh
source ~/projects/scripts/bash/git_prompt.sh

# old PS1
# export PS1="\h:\W \u\$ "
export PS1="\[\033[36m\]\u\[\033[m\]@\[\033[32m\]\h:\[\033[33;1m\]\w\[\033[m\]\$(__git_ps1 '(%s)')$ "
export CLICOLOR=1
export LSCOLORS=ExFxBxDxCxegedabagacad
alias ls='ls -GFh'

if [ -f ~/.bash_variables ]; then
  . ~/.bash_variables
fi

if [ -f ~/.bash_functions ]; then
  . ~/.bash_functions
fi

export PATH=${PATH}:${USER_SCRIPTS}:${HOME}/.local/bin