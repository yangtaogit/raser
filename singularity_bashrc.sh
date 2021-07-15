# BashRC file for RASER inside Singularity aser 
# Author SHI Xin <shixin@ihep.ac.cn>
# Created [2021-06-30 Wed 14:25]


#--------------------------------------------------
# Prompt 
#--------------------------------------------------
HOST=$(hostname | cut -d. -f1)
if [ "$PS1" ]; then
  PS1="\u_$HOST RASER> "
  ignoreeof=1
fi


#--------------------------------------------------
# ENV settings
#--------------------------------------------------
source $HOME/geant4/install/bin/geant4.sh 
export PYTHONPATH=$PYTHONPATH:$HOME/geant4/install/lib/python3.8/site-packages
export GEANT4_INSTALL=$HOME/geant4/install
export DISPLAY=localhost:11.0

#--------------------------------------------------
# Aliases
#--------------------------------------------------
alias l="ls --color"
alias ll="l -lh"
alias p="pwd"
alias rl="root -l" 

