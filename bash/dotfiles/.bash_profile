###########
###  START Custom Profile
# Set proxies for use on the terminal
if [ `ifconfig | grep -c "141.221.120\|141.221.121\|141.221.124\|141.221.113"` -ge 1 ]; then
  echo 'Using FN Proxies'
  proxy="http://fnwebbalance.fn.inl.gov:8080"
#elif [ `ifconfig | grep -c "134.20.\|141.221."` -ge 1 ]; then
#  echo 'Using INL Proxies'
#  proxy="http://webbalance.inl.gov:8080"
else
  # No Proxies
  echo 'No Proxies set'
  proxy=""
fi

if [ "x$proxy" != "x" ]; then
  proxy_array=(http_proxy https_proxy ftp_proxy ftps_proxy)
  for proxy_item in ${proxy_array[*]}; do
    export $proxy_item=$proxy
  done
fi
# End Set proxies
###

###
# Set Pretty stuff
#PS1="\[\033[1;34m\][\u]\[\033[1;32m\][\w]\[\033[0m\]> "
#export CLICOLOR=1
#export LSCOLORS=exfxcxdxbxegedabagacad
#
###

###
# Set useful aliases
alias ll="ls -latrh"
alias l="ls -latrh"
# End aliases
###
###  END Custom Profile
###  you can safely append any custom content below this line

# For autojump capabilities
[ -f /usr/local/etc/profile.d/autojump.sh ] && . /usr/local/etc/profile.d/autojump.sh

# Custom settings from Jarin - stored in .bashrc
if [ -f ~/.bashrc ]; then
  . ~/.bashrc
fi

#MOOSE_ENVIRONMENT
# Uncomment to enable MOOSE prompt:
# export MOOSE_PROMPT=true

# Uncomment to enable autojump:
# export MOOSE_JUMP=true

# Source MOOSE profile
if [ -f /opt/moose/environments/moose_profile ]; then
        . /opt/moose/environments/moose_profile

        # Make the moose compiler stack available.
        # Note: if you have any additional package managers installed
        # (Homebrew, Macports, Fink, etc) know that you must perform
        # the following command _after_ the source commands for the
        # afore mentioned package managers.
        module load moose-dev-clang
fi
#ENDMOOSE_ENVIRONMENT
