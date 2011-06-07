#!/bin/zsh 


if [ "$BASEDIR" = "" ]; then 
  export BASEDIR=$PWD
fi

echo BASEDIR $BASEDIR

# How to access the svn directory
# svn co svn+ssh://login.hepforge.org/hepforge/svn/applgrid/
# svn co http://svn.hepforge.org/applgrid/trunk applgrid
# chmod u+x appl_grid/bin/applgrid-config
#
#
# SLC5 set-up 64-bit machine
#
#export platform=x86_64-slc5-gcc43-opt
#export ROOTSYS=/afs/cern.ch/sw/lcg/app/releases/ROOT/5.26.00b/slc4_amd64_gcc34/root
#export LHAPDFPATH=/afs/cern.ch/sw/lcg/external/MCGenerators/lhapdf/5.8.3/slc4_amd64_gcc34
#export FASTJETPATH=/afs/cern.ch/sw/lcg/external/fastjet/2.4.2/slc4_amd64_gcc34
#
#
#export PATH=${BASEDIR}/bin:${ROOTSYS}/bin:${LHAPDFPATH}/bin:${FASTJETPATH}/bin:${PATH}
#export LD_LIBRARY_PATH=/usr/lib/:${BASEDIR}/lib:${ROOTSYS}/lib:${LHAPDFPATH}/lib:${FASTJETPATH}/lib
#
#
export platform=x86_64-slc5-gcc43-opt
#
#export external=/afs/cern.ch/sw/lcg/external 
#export releases=/afs/cern.ch/sw/lcg/app/releases
export contrib=/afs/cern.ch/sw/lcg/contrib
#
export GCC43_INSTALL_DIR=${contrib}/gcc/4.3.2/${platform}
source ${contrib}/gcc/4.3/x86_64-slc5/setup.sh
export PATH=${GCC43_INSTALL_DIR}/bin:${PATH}
export LD_LIBRARY_PATH=${GCC43_INSTALL_DIR}/lib64:${LD_LIBRARY_PATH}
#
export CLHEPDIR=/afs/cern.ch/sw/lcg/external/clhep/2.0.4.0/${platform}
export ROOTSYS=/afs/cern.ch/sw/lcg/app/releases/ROOT/5.28.00b/${platform}/root
export LHAPDFPATH=/afs/cern.ch/sw/lcg/external/MCGenerators/lhapdf/5.8.5/${platform}

##
export FASTJETPATH=/afs/cern.ch/sw/lcg/external/fastjet/2.4.2/${platform}
#
##
##  longer AFS token
##
##  reauth 1000000 pavel --fork
#
export LD_LIBRARY_PATH=$HOME/lib:$CLHEPDIR/lib:$ROOTSYS/lib:$LHAPDFPATH/lib/:$FASTJETPATH/lib/:${HOME}/lib/EoP:$LD_LIBRARY_PATH

export PATH=$ROOTSYS/bin:$CLHEPDIR/bin:$HOME/bin:$FASTJETPATH/bin:$LHAPDFPATH/bin:$PATH
##
export BASEDIR=$PWD
export ARCH=`$BASEDIR/bin/rootarch.sh`
export CXXFLAGS=$ARCH
export F90FLAGS=$ARCH
export CFLAGS=$ARCH
export FFLAGS=$ARCH
export LDFLAGS=$ARCH

export PATH=${BASEDIR}/bin:${PATH}
export LD_LIBRARY_PATH=${BASEDIR}/lib:${LD_LIBRARY_PATH}


##   Fortran 95
##
#LHDIR=/afs/cern.ch/sw/fortran/lahey/lf9562c
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$LHDIR/lib
#export PATH=$LHDIR/bin:$PATH
#export MANPATH=${MANPATH}:$LHDIR/manuals/man/lf95



