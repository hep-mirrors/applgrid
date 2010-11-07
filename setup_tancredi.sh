#!/bin/zsh 


if [ "$BASEDIR" = "" ]; then 
  export BASEDIR=$PWD
fi


# How to access the svn directory
#svn co svn+ssh://login.hepforge.org/hepforge/svn/applgrid/
# svn co http://svn.hepforge.org/applgrid/trunk applgrid
# chmod u+x appl_grid/bin/applgrid-config
#
# Ubunuty 9 set-up on my machine
#export ROOTSYS=/home/tcarli/root/root5.24
#export LHAPDF=/home/tcarli/lhapdf-5.7.0
#
#
# SLC4 set-up, 64-bit machine
#export ROOTSYS=/afs/cern.ch/sw/lcg/app/releases/ROOT/5.24.00/slc4_amd64_gcc34/root
#export LHAPDF=/afs/cern.ch/sw/lcg/external/MCGenerators/lhapdf/5.7.1/slc4_amd64_gcc34
#
#
# SLC5 set-up 43-bit machine
#
#export ROOTSYS=/afs/cern.ch/sw/lcg/app/releases/ROOT/5.24.00/x86_64-slc5-gcc34-dbg/root
#export FASTJETPATH=/afs/cern.ch/sw/lcg/external/fastjet/2.4.2/slc4_amd64_gcc34
#
#export ROOTSYS=/afs/cern.ch/sw/lcg/app/releases/ROOT/5.24.00/x86_64-slc5-gcc43-dbg/root
#export LHAPDF=/afs/cern.ch/sw/lcg/external/MCGenerators/lhapdf/5.8.4/x86_64-slc5-gcc43-opt
#export FASTJETPATH=/afs/cern.ch/sw/lcg/external/fastjet/2.4.2/x86_64-slc5-gcc43-opt

export ROOTSYS=/afs/cern.ch/sw/lcg/app/releases/ROOT/5.26.00b/slc4_amd64_gcc34/root
export LHAPDFPATH=/afs/cern.ch/sw/lcg/external/MCGenerators/lhapdf/5.8.3/slc4_amd64_gcc34
export FASTJETPATH=/afs/cern.ch/sw/lcg/external/fastjet/2.4.2/slc4_amd64_gcc34
#
#
#export PATH=/usr/sue/bin:/usr/local/bin:/usr/local/bin/X11:/usr/bin:/bin:/usr/bin/X11
export PATH=${BASEDIR}/bin:${ROOTSYS}/bin:${LHAPDFPATH}/bin:${FASTJETPATH}/bin:${PATH}
export LD_LIBRARY_PATH=/usr/lib/${BASEDIR}/lib:${LHAPDFPATH}/lib:${FASTJETPATH}/lib
#
#
#
# this adds lhapdf/bin to path so that the essential  
# lhapdf-config script can be found
export PATH=$ROOTSYS/bin:$LHAPDF/bin:$PATH
export LD_LIBRARY_PATH=$ROOTSYS/lib:$LD_LIBRARY_PATH
#
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$LHAPDF/lib
export PATH=$PATH:/cern/pro/bin:$LHAPDF/bin:$BASEDIR/appl_grid/bin:$BASEDIR/bin
#
#
# LHDIR=/afs/cern.ch/sw/fortran/lahey/lf9562c
# export LD_LIBRARY_PATH=$LHDIR/lib
# export PATH=$LHDIR/bin:$PATH
# export MANPATH=${MANPATH}:$LHDIR/manuals/man/lf95

