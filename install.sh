#!/bin/sh -v

if [ "$COPYDIR" = "" ]; then 
  export COPYDIR=$PWD
fi

echo "copy directory " $COPYDIR


if [ "$BASEDIR" = "" ]; then 
  export BASEDIR=$PWD
fi
export INSTALLBASE=$BASEDIR

echo "base directory " $BASEDIR 
if [ -e $BASEDIR ]; then 
   echo "$BASEDIR already exists"
else
   mkdir $BASEDIR
fi

cd $BASEDIR


pwd


export ARCH=`rootarch.sh`

export CXXFLAGS=$ARCH
export F90FLAGS=$ARCH
export F77FLAGS=$ARCH
export CFLAGS=$ARCH
export FFLAGS=$ARCH
export LDFLAGS=$ARCH



#tar -xzf $COPYDIR/appl_grid.tgz
#tar -xzf $COPYDIR/pdf-conv-1.0.tgz
## tar -xzf $COPYDIR/mcfm.tgz mcfm
#tar -xzf $COPYDIR/lhpdf-1.0.0.tgz 
#tar -xzf $COPYDIR/nlojet++-4.0.1.tgz 
#tar -xzf $COPYDIR/ap.tgz	
#
# exit

export PATH=$BASEDIR/bin:$PATH
export LD_LIBRARY_PATH=$BASEDIR/lib:$BASEDIR/libexec:$LD_LIBRARY_PATH
export DYLD_LIBRARY_PATH=$BASEDIR/lib:$BASEDIR/libexec:$DYLD_LIBRARY_PATH



############################
# INSTALL THE APPL GRID CODE 
# AND PDF EVOLUTION CODE
############################
cd $BASEDIR/appl_grid/src

if [ "$1" = "clean" ]; then
 make clean
fi


make install



cd $BASEDIR/pdf-conv-1.0/src

if [ "$1" = "clean" ]; then
 make clean
fi

make  install

######################
# INSTALL MCFM
######################

cd $BASEDIR

# if [ ! -e $BASEDIR/lib/libgfortran.so ]; then
#    ln -s /usr/lib/libgfortran.so.2 $BASEDIR/lib/libgfortran.so
# fi
#

############################################
# BUILD MCFM TEST EXECUTABLE
############################################
cd $BASEDIR/mcfm
make install
#
echo "Don't forget to set the LD_LIBRARY_PATH for root and lhapdf if you have to"
#

############################################
# RUN THE MCFM TEST EXECUTABLE
############################################

cd $BASEDIR/mcfm/run
#
../exe/*/mcfm Winput.DAT >&  mcfm-0.log
../exe/*/mcfm Winput.DAT >&  mcfm-1.log
../exe/*/stand grid-30-Wweight_eta4.root
#
# echo EXITING......
# exit


#############################
# INSTALL NLOJET PDF LIBRARY 
#############################

cd $BASEDIR


echo  "bld-pdf" 

if [ "$1" = "clean" ]; then
    rm -rf bld-pdf
fi

if [ "$1" = "nloclean" ]; then
    rm -rf bld-pdf
fi


if [ -e bld-pdf ]; then 
  echo "bld-pdf already exists" 
else

  mkdir bld-pdf

  if [ "$1" = "clean" ]; then
    make clean
  fi

  cd bld-pdf
  ../lhpdf-1.0.0/configure --prefix=$BASEDIR
  # > $BASEDIR/lhpdf.log

  make        
  #  >> $BASEDIR/lhpdf.log

  make install 
  #  >> $BASEDIR/lhpdf.log

fi


######################
# INSTALL NLOJET
######################

cd $BASEDIR


echo "bld-nlo" 
if [ -e bld-nlo ]; then 
  echo "bld-nlo already exists" 
else
  mkdir bld-nlo

  cd   bld-nlo

  echo "nlojet"
  ../nlojet++-4.0.1/configure --prefix=$BASEDIR 
  # cd blo> $BASEDIR/nlo.log

  if [ "$1" = "clean" ]; then
    make clean
  fi

  make          
  # >> $BASEDIR/nlo.log
  make install 
  # >> $BASEDIR/nlo.log

fi

echo $PATH
echo "which create-nlojet-user"
#echo "why. oh why, do we need these 'helper' routines that make automating things so difficult..."
which create-nlojet-user



######################
# THREEJET MODULES
######################

# cd $BASEDIR

# echo "modules"
# if [ -e modules-4.0.0.tar.gz  ]; then 
#   echo "modules-4.0.0 already exists"
# else
#   cp $COPYDIR/modules-4.0.0.tar.gz modules-4.0.0.tar.gz
# fi

# tar -zxf modules-4.0.0.tar.gz

# cd "modules-4.0.0/pp->jets/hep-ph:0307268/ThreeJet"
# make all  


#  nlojet++ --calculate -u kT-scale.la
#  nlojet++ --calculate -u kT-xsec.la
#  nlojet++ --calculate -u cone-scale.la


######################
# INSTALL GRID FILLING
# MODULES 
######################

cd $BASEDIR


echo "appl-paper"
if [ -e $BASEDIR/appl-paper  ]; then
   echo "appl-paper already exists"
else 
#   cp $COPYDIR/ap.tgz ap.tgz
    echo doing nowt
fi

	
echo "building fillgrid..."
cd appl-paper

if [ "$1" = "clean" ]; then
 make clean
fi

make fillgrid

make standalone

echo PATH is: $PATH

echo LD_LIBRARY_PATH: $LD_LIBRARY_PATH

if [ -e PDFsets ]; then 
  ls -ld PDFsets
else
  ln -s `lhapdf-config --pdfsets-path` .
  ls -ld PDFsets
fi



which nlojet++
# nlojet++ --calculate -u fillgrid.la --max-event 20000
nlojet++ --calculate -u fillgrid.la 
