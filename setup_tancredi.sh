#!/bin/zsh 

# How to access the svn directory
export SVNUSR=svn+ssh://svn.cern.ch/reps/atlasusr/tcarli

#export ROOTSYS=/home/tcarli/root/root5.24
export LHAPDF=/home/tcarli/lhapdf-5.7.0


export PATH=/usr/sue/bin:/usr/local/bin:/usr/local/bin/X11:/usr/bin:/bin:/usr/bin/X11
export LD_LIBRARY_PATH=/usr/lib

export PATH=$ROOTSYS/bin:$LHAPDF/bin:$PATH
export LD_LIBRARY_PATH=$ROOTSYS/lib:$LD_LIBRARY_PATH

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$LHAPDF/lib
export PATH=$PATH:/cern/pro/bin:$LHAPDF/bin:$BASEDIR/bin

