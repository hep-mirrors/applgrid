#!/bin/sh 

######################################################
#  installation script - now greatly improved 
#  by default installs and runs everything 
#  except mcfm
#  
#  has control flags to allow any combination 
#  of jobs to be done. 
# 
#    APPLGRID - install appl_grid
#    PDF      - install hoppet
#    MCFM     - install mcfm
#    RMCFM    - run mcfm
#    NLO      - install nlojet
#    NLOMOD   - instal nlojet module
#    RNLOMOD  - run nlojet module
#  
#  To do a single job, 
#  diasble everything with --none, ie 
#
#   install.sh --none --appl
#
######################################################


ARGS=$*

# control flags

APPLGRID=1
PDF=1
MCFM=0
RMCFM=0
NLO=1
NLOMOD=1
RNLOMOD=1


#   echo "setting everything"
setall() { 
    APPLGRID=1
    PDF=1
    MCFM=1
    RMCFM=1
    NLO=1
    NLOMOD=1
    RNLOMOD=1
}

#   echo "unsetting everything"
unsetall() { 
    APPLGRID=0
    PDF=0
    MCFM=0
    RMCFM=0
    NLO=0
    NLOMOD=0
    RNLOMOD=0
}

# usage message
usage() { 
    echo "usage: install.sh [OPTIONS]\n"
    echo "  --help| -h   this help"
    echo "  --all|-a     do everything"
    echo "  --none|-n    don't do anything"
    echo "  --appl       install appl_grid"
    echo "  --pdf        install hoppet"
    echo "  --mcfm       install mcfm"
    echo "  --nlo        install nlojet"
    echo "  --mod        install nlojet module"
    echo "  --runmcfm    run mcfm"
    echo "  --runmod     run nlojet module"
    echo "\nReport bugs to <sutt@cern.ch>" 
#   exit
}


###############################################
# parse arguments and set control flags 
###############################################

for WORD in $ARGS ; do
   case $WORD in
       --all|-a)   setall;;  
       --none|-n)  unsetall;;
       --appl)     APPLGRID=1;;
       --pdf)      PDF=1;;
       --mcfm)     MCFM=1;;
       --runmcfm)  RMCFM=1;;
       --nlo)      NLO=1;;
       --mod)      NLOMOD=1;;
       --runmod)   RNLOMOD=1;;
       --help|-h)  usage;exit;;
   esac
done



###############################################
#  setup everything
###############################################

if [ "$BASEDIR" = "" ]; then 
  export BASEDIR=$PWD
fi
export INSTALLBASE=$BASEDIR

echo "base directory " $BASEDIR 
if [ -e $BASEDIR ]; then 
   echo "installing in $BASEDIR"
else
   echo "$BASEDIR does not exists" 
   exit
fi

cd $BASEDIR


export PATH=$BASEDIR/bin:$PATH

ROOTFOUND=`which root-config`

if [ "$ROOTFOUND" = "" ]; then 
  echo "root not setup"
  exit
fi

LHAPDFFOUND=`which lhapdf-config`
 
if [ "$LHAPDFFOUND" = "" ]; then 
  echo "lhapdf not setup"
  exit
fi



#  use root-config to decide what architecture root libraries will be used
 
export ARCH=`rootarch.sh`

export CXXFLAGS=$ARCH
export F90FLAGS=$ARCH
export F77FLAGS=$ARCH
export CFLAGS=$ARCH
export FC=gfortran
# export FFLAGS="$ARCH -O"
export FFLAGS=$ARCH
export LDFLAGS=$ARCH




export ROOTSYS=`root-config --prefix`
export LHAPDFLIB=`lhapdf-config --prefix`/lib


export PATH=$ROOTSYS/bin:$PATH
export LD_LIBRARY_PATH=$ROOTSYS/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$BASEDIR/lib:$BASEDIR/libexec:$LHAPDFLIB:$LD_LIBRARY_PATH
export DYLD_LIBRARY_PATH=$BASEDIR/lib:$BASEDIR/libexec:$LHAPDFLIB:$DYLD_LIBRARY_PATH




############################
# INSTALL THE APPL GRID CODE 
############################
install_appl_grid() { 

    cd $BASEDIR/appl_grid/src

    if [ "$1" = "clean" ]; then
	make clean
    fi

    make install

}



############################
# INSTALL EVOLUTION CODE
############################
install_pdf() { 

    cd $BASEDIR/hoppet
    if [ "$1" = "clean" ]; then
	make clean
    fi
    ./configure --prefix=$BASEDIR FC=$FC FFLAGS="$FFLAGS" LDFLAGS="$LDFLAGS"
    make
    make check
    make install

    # cd $BASEDIR/pdf-conv-1.0/src
    # 
    # if [ "$1" = "clean" ]; then
    # 	make clean
    # fi
    # 
    # make  install
    
}



######################
# INSTALL MCFM
######################
install_mcfm() { 
    cd $BASEDIR
    cd $BASEDIR/mcfm
    make install    
}


############################################
# RUN THE MCFM TEST EXECUTABLE
############################################
run_mcfm() { 

    cd $BASEDIR/mcfm/run
    rm -f *.log

    if [ -e ../exe/$ARCH/mcfm ]; then 
      ../exe/$ARCH/mcfm Winput.DAT >  mcfm-0.log
      ../exe/$ARCH/mcfm Winput.DAT >  mcfm-1.log
      ../exe/$ARCH/stand grid-30-Wweight_eta4.root
    fi
}


install_nlojet() { 
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
	
	make        
	make install 
	
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
	
	if [ "$1" = "clean" ]; then
	    make clean
	fi
	
	make          
	make install 
	
    fi
    
    echo $PATH
    echo "which create-nlojet-user"

#   which create-nlojet-user
    
}




######################################
# INSTALL NLOJET GRID FILLING MODULES 
######################################
install_nlojet_module() { 

    echo "jetmod"
    if [ -e $BASEDIR/jetmod  ]; then
	echo "cd jetmod"
	cd $BASEDIR/jetmod
    else 
	echo "$BASEDIR/jetmod doesn't exist"
	exit
    fi
    

    echo "building fillgrid..."
    
    if [ "$1" = "clean" ]; then
	make clean
    fi
    
    make install
    
    echo PATH is: $PATH
    echo LD_LIBRARY_PATH: $LD_LIBRARY_PATH
    
    if [ -e PDFsets ]; then 
	ls -ld PDFsets
    else
	ln -s `lhapdf-config --pdfsets-path` .
	ls -ld PDFsets
    fi
    
}


##################################
# RUN NLOJET GRID FILLING MODULES 
##################################
run_nlojet_module () { 

    cd $BASEDIR/jetmod

    if [ -e output/weight_c.root ]; then
	echo "nlojet already executed"
    else  
	which nlojet++
#       nlojet++ --calculate -u fillgrid.la --max-event 20000
#       generate the grid phase space
        nlojet++ --calculate -c full -u fillgrid.la --max-event 1000010
#       now optimise it and generate the weights properly
        nlojet++ --calculate -c full -u fillgrid.la --max-event 1000010
    fi

    if [ -e output/weight_c.root ]; then 
#	./standalone output/weight_c.root
	./stand      output/weight_c.root
    else
	echo "nlojet grid code did not complete correctly" 
    fi
    
}


###########################
# ACTUALLY RUN THE STAGES 
###########################

if [ "$PDF" = 1 ];      then install_pdf        $ARGS; fi
if [ "$APPLGRID" = 1 ]; then install_appl_grid  $ARGS; fi
if [ "$MCFM" = 1 ];     then install_mcfm       $ARGS; fi
if [ "$RMCFM" = 1 ];    then run_mcfm;                 fi        
if [ "$NLO" = 1 ];      then    install_nlojet         $ARGS; fi
if [ "$NLOMOD" = 1 ];   then    install_nlojet_module  $ARGS; fi
if [ "$RNLOMOD" = 1 ];  then    run_nlojet_module;            fi 
