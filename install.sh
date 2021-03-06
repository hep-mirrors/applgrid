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
#    MCFM60   - install mcfm 6.0
#    RMCFM    - run mcfm
#    RMCFM60  - run mcfm 6.0
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


###############################################
#  setup everything
###############################################

if [ "$BASEDIR" = "" ]; then 
  export BASEDIR=$PWD
fi
export INSTALLBASE=$BASEDIR


# control flags

APPLGRID=1
PDF=1
MCFM=0
MCFM60=1
RMCFM=0
RMCFM60=0
NLO=1
NLOMOD=1
RNLOMOD=1
FASTJET=1
USERJ=1

#   echo "setting everything"
setall() { 
    APPLGRID=1
    PDF=1
    MCFM=1
    MCFM60=1
    RMCFM=1
    RMCFM60=1
    NLO=1
    NLOMOD=1
    RNLOMOD=1
    FASTJET=1
    USERJ=1
}

#   echo "unsetting everything"
unsetall() { 
    APPLGRID=0
    PDF=0
    MCFM=0
    MCFM60=0
    RMCFM=0
    RMCFM60=0
    NLO=0
    NLOMOD=0
    RNLOMOD=0
    FASTJET=0
    USERJ=0
}


#   echo "setting everything"
printflags() { 
    echo "APPLGRID $APPLGRID"
    echo "PDF     $PDF"
    echo "MCFM    $MCFM"
    echo "MCFM60  $MCFM60"
    echo "RMCFM   $RMCFM"
    echo "RMCFM60 $RMCFM60"
    echo "NLO     $NLO"
    echo "NLOMOD  $NLOMOD"
    echo "RNLOMOD $RNLOMOD"
    echo "FASTJET $FASTJET"
    echo "USERJ   $USERJ"
}


# usage message
usage() { 
    echo "usage: install.sh [OPTIONS]\n"
    echo "  --help| -h   this help"
    echo 
    echo "  --all|-a        do everything"
    echo "  --none|-n       don't do anything"
    echo "  --appl          install appl_grid"
    echo "  --pdf|--hoppet  install hoppet"
    echo "  --mcfm          install mcfm"
    echo "  --mcfm60        install mcfm v6.0"
    echo "  --nlo           install nlojet"
    echo "  --mod           install nlojet module"
    echo
    echo "  --runmcfm       run mcfm"
    echo "  --runmcfm60     run mcfm v6.0"
    echo "  --runmod        run nlojet module"
    echo 
    echo "  --user          build and run simple user example"
    echo "  --fastjet       build fastjet"
    echo
    echo "  --prefix=<directory>  set the installation directory"
    echo 
    echo "\nReport bugs to <sutt@cern.ch>" 
#   exit
}


# set the installation diectory if required
prefix() {
    export INSTALLBASE=`echo $PREFIX | sed "s#--prefix=##"`
#   echo installing in $INSTALLBASE
} 



###############################################
# parse arguments and set control flags 
###############################################


for WORD in $ARGS ; do
   case $WORD in
       --all|-a)   setall;;  
       --none|-n)  unsetall;;
       --appl)     APPLGRID=1;;
       --pdf|--hoppet)      PDF=1;;
       --mcfm)     MCFM=1;;
       --mcfm60)   MCFM60=1;;
       --runmcfm)  RMCFM=1;;
       --runmcfm60) RMCFM60=1;;
       --nlo)      NLO=1;;
       --mod)      NLOMOD=1;;
       --runmod)   RNLOMOD=1;;
       --user)     USERJ=1;;
       --fastjet)  FASTJET=1;;
       --prefix=*) PREFIX=$WORD;prefix;;
       --help|-h)  usage;exit;;
       *)          usage;exit;;
   esac
done


echo "base directory " $BASEDIR 
if [ -e $INSTALLBASE ]; then 
   echo "installing in $INSTALLBASE"
else
   echo "$INSTALLBASE does not exist" 
   exit
fi


cd   $BASEDIR

export PATH=$INSTALLBASE/bin:$BASEDIR/bin:$PATH

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
 
# export ARCH=`rootarch.sh`
export ARCH=`root-config --ldflags`
export SYSARCH=`setarch.sh`

export CXXFLAGS=$ARCH
export F90FLAGS=$ARCH
export F77FLAGS=$ARCH
export CFLAGS=$ARCH
export FC=gfortran
# export FFLAGS="$ARCH -O"
export FFLAGS=$ARCH
export LDFLAGS=$ARCH

export HOPPETLIBS=
export HOPPETINCS=
export HOPPETFLAG=  

export FASTJETLIBS=
export FASTJETINCS=
export FASTJETFLAG=  




export ROOTSYS=`root-config --prefix`
export LHAPDFLIB=`lhapdf-config --prefix`/lib


export PATH=$ROOTSYS/bin:$PATH
export LD_LIBRARY_PATH=$ROOTSYS/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$INSTALLBASE/lib:$INSTALLBASE/libexec:$LHAPDFLIB:$LD_LIBRARY_PATH
export DYLD_LIBRARY_PATH=$INSTALLBASE/lib:$INSTALLBASE/libexec:$LHAPDFLIB:$DYLD_LIBRARY_PATH




############################
# INSTALL THE APPL GRID CODE 
############################
install_appl_grid() { 

    cd $BASEDIR/appl_grid

    ./configure --prefix=$INSTALLBASE

    #    if [ "$1" = "clean" ]; then
    #	   make clean
    #    fi

    make clean
    autoreconf -i
 
    make 
    make install

    hash -r

}



############################
# INSTALL EVOLUTION CODE
############################
install_pdf() { 

    cd $BASEDIR/hoppet

    if [ -e $BASEDIR/hoppet/benchmarking/benchmark.f90 ]; then 
         rm $BASEDIR/hoppet/benchmarking/benchmark.f90 
    fi

    if [ "$1" = "clean" ]; then
	make clean
    fi
    ./configure --prefix=$INSTALLBASE FC=$FC FFLAGS=" -fPIC $FFLAGS" LDFLAGS="$LDFLAGS"
    make
#   make check
    make install

    hash -r

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

install_mcfm60() {
    cd $BASEDIR
    cd $BASEDIR/mcfm-6.0
    make install    
#   should no longer be needed by lhapdf. hmmmm
    cd Bin
    if [ -e PDFsets ]; then 
	ls -ld PDFsets
    else
	ln -s `lhapdf-config --pdfsets-path` .
	ls -ld PDFsets
    fi
    

}


############################################
# RUN THE MCFM TEST EXECUTABLE
###########################cd ..#################
run_mcfm() { 

  if [ -d $BASEDIR/mcfm/run ]; then 

    cd $BASEDIR/mcfm/run
    rm -f *.log

#   should no longer be needed by lhapdf. hmmmm
    if [ -e PDFsets ]; then 
	ls -ld PDFsets
    else
	ln -s `lhapdf-config --pdfsets-path` .
	ls -ld PDFsets
    fi
    
    echo "running mcfm - please wait ..."

    if [ -e ../exe/$SYSARCH/mcfm ]; then 

       rm -rf W_*	

       #### Wminus test executable ####
#      ../exe/$SYSARCH/mcfm WMinput.DAT >&  mcfm-wm0.log
#      ../exe/$SYSARCH/mcfm WMinput.DAT >&  mcfm-wm1.log
#      ../exe/$SYSARCH/stand grid-30-Wminus_eta4.root
#      mv fout.root fout-Wminus.root

      rm -rf W_*	

#     #### Wplus test executable ####
      ../exe/$SYSARCH/mcfm WPinput.DAT >&  mcfm-wp0.log
      ../exe/$SYSARCH/mcfm WPinput.DAT >&  mcfm-wp1.log
      ../exe/$SYSARCH/stand grid-30-Wplus_eta4.root
      mv fout.root fout-Wplus.root

      cd ..
    fi
  fi
}


run_mcfm60() { 

  if [ -d $BASEDIR/mcfm-6.0/Bin ]; then 

    cd $BASEDIR/mcfm-6.0/Bin
    rm -f *.log

#   should no longer be needed by lhapdf. hmmmm
    if [ -e PDFsets ]; then 
	ls -ld PDFsets
    else
	ln -s `lhapdf-config --pdfsets-path` .
	ls -ld PDFsets
    fi
    
    echo "running mcfm - please wait ..."

    if [ -e mcfm ]; then 

      rm -rf W_*	

#     #### Wplus test executable ####
      mcfm winput.DAT >&  mcfm.log
      mcfm winput.DAT >&  mcfm.log
      stand grid-30-Wplus_eta4.root
      mv fout.root fout-Wplus.root

      cd ..
    fi
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
	../lhpdf-1.0.0/configure --prefix=$INSTALLBASE
	
	make        
	make install 

    fi

    hash -r


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
	../nlojet++-4.0.1/configure --prefix=$INSTALLBASE 
	
	if [ "$1" = "clean" ]; then
	    make clean
	fi
	
	make          
	make install 
	
    fi
    
    hash -r

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
    
    hash -r

    echo PATH is: $PATH
    echo LD_LIBRARY_PATH: $LD_LIBRARY_PATH
 
    # recent LHAPDF should no longer need this   
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


run_user () { 
    cd $BASEDIR/user

    if [ "$1" = "clean" ]; then
	make clean
    fi

    make all

    echo "\nrunning executables..."

    if [ -e ../jetmod/output/weight_c.root ]; then 
      ./stand ../jetmod/output/weight_c.root
    else
       echo "weight file not found (../jetmod/output/weight_c.root)"
       echo "have you already run nlojet++ ?"
    fi
       
    if [ -e fnl0004.tab ]; then 
      ./fnlo  fnl0004.tab
    else 
      echo "fastnlo grid file (fnl0004.tab) missing" 	
    fi

}

install_fastjet () { 

    cd $BASEDIR

    if [ -e "$BASEDIR/fastjet-*" ]; then
       echo "fastjet already installed"
    else 
      FASTJETVER=`wget --quiet http://www.lpthe.jussieu.fr/~salam/fastjet/version.html -O - | grep '^[0-9]' -` 
      curl http://www.lpthe.jussieu.fr/~salam/fastjet/repo/fastjet-$FASTJETVER.tar.gz > /tmp/fastjet-$FASTJETVER.tar.gz
      tar -xzf /tmp/fastjet-$FASTJETVER.tar.gz 

      if [ -e "$BASEDIR/fastjet-$FASTJETVER" ]; then 
         cd $BASEDIR/fastjet-$FASTJETVER
         ./configure --prefix=$INSTALLBASE  FC=$FC CXXFLAGS="$CXXFLAGS" FFLAGS="$FFLAGS" CFLAGS="$CFLAGS" LDFLAGS="$LDFLAGS"
         if [ "$1" = "clean" ]; then
 	   make clean
         fi
         make install
      fi

    fi

}


###########################
# ACTUALLY RUN THE STAGES 
###########################

hash -r

printflags

if [ "$PDF" = 1 ];      then install_pdf        $ARGS; fi

if [ -e $INSTALLBASE/bin/hoppet-config ]; then 
    export HOPPETLIBS="` $INSTALLBASE/bin/hoppet-config --libs `"
    export HOPPETINCS="` $INSTALLBASE/bin/hoppet-config --cxxflags `" 
    export HOPPETFLAG="  -DHOPPET " 
fi 

if [ "$APPLGRID" = 1 ]; then install_appl_grid  $ARGS; fi
if [ "$MCFM" = 1 ];     then install_mcfm       $ARGS; fi
if [ "$MCFM60" = 1 ];   then install_mcfm60     $ARGS; fi
if [ "$RMCFM" = 1 ];    then run_mcfm;                 fi        
if [ "$RMCFM60" = 1 ];  then run_mcfm60;               fi        
if [ "$NLO" = 1 ];      then install_nlojet     $ARGS; fi


FASTJETFOUND=`which fastjet-config`

if [ "$FASTJETFOUND" = "" ]; then  
  if [ "$FASTJET" = 1 ];  then    install_fastjet        $ARGS; fi
 
  if [ -e $INSTALLBASE/bin/fastjet-config ]; then 
     export FASTJETLIBS="` $INSTALLBASE/bin/fastjet-config --libs `"
     export FASTJETINCS="` $INSTALLBASE/bin/fastjet-config --cxxflags `" 
     export FASTJETFLAG="  -DFASTJET " 
  fi
else
  if [ -e $INSTALLBASE/bin/fastjet-config ]; then 
     export FASTJETLIBS="` fastjet-config --libs `"
     export FASTJETINCS="` fastjet-config --cxxflags `" 
     export FASTJETFLAG="  -DFASTJET " 
  fi
fi 

if [ "$NLOMOD" = 1 ];   then    install_nlojet_module  $ARGS; fi
if [ "$RNLOMOD" = 1 ];  then    run_nlojet_module;            fi 
if [ "$USERJ" = 1 ];    then    run_user           $ARGS;  fi 
