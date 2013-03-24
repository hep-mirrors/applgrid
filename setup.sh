
#
#-----------------------
#  setup for SLC5 
#
#export platform=x86_64-slc5-gcc43-opt
#
#export external=/afs/cern.ch/sw/lcg/external
#export releases=/afs/cern.ch/sw/lcg/app/releases
#export contrib=/afs/cern.ch/sw/lcg/contrib
#
#
#export GCC43_INSTALL_DIR=${contrib}/gcc/4.3.2/${platform}
#source ${contrib}/gcc/4.3/x86_64-slc5/setup.sh
#export PATH=${GCC43_INSTALL_DIR}/bin:${PATH}
#export LD_LIBRARY_PATH=${GCC43_INSTALL_DIR}/lib64:${LD_LIBRARY_PATH}
#
#export ROOTSYS=/afs/cern.ch/sw/lcg/app/releases/ROOT/5.28.00b/x86_64-slc5-gcc43-opt/root
#export PATH=${PATH}:${ROOTSYS}/bin
#export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${ROOTSYS}/lib
#
#------------------------------
#                                                                                                                                


export BASEDIR=$PWD
export ARCH=`$BASEDIR/bin/rootarch.sh`
export CXXFLAGS=$ARCH
export F90FLAGS=$ARCH
export CFLAGS=$ARCH
export FFLAGS=$ARCH
export LDFLAGS=$ARCH

export PATH=${BASEDIR}/bin:${PATH}
export LD_LIBRARY_PATH=${BASEDIR}/lib:${LD_LIBRARY_PATH}

# to redo the makefiles etc
#autoreconf --install --force

