#!/bin/sh -v


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
export FFLAGS="$ARCH -O"
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
    ./configure --prefix=$BASEDIR FC=$FC FFLAGS="$FFLAGS"
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
    rm *.log

    ../exe/*/mcfm Winput.DAT >  mcfm-0.log
    ../exe/*/mcfm Winput.DAT >  mcfm-1.log
    ../exe/*/stand grid-30-Wweight_eta4.root
    
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




######################
# INSTALL GRID FILLING
# MODULES 
######################
install_nlojet_module() { 

    echo "nlojet-module"
    if [ -e $BASEDIR/nlojet-module  ]; then
	echo "cd nlojet-module"
	cd $BASEDIR/nlojet-module
    else 
	echo "$BASEDIR/nlojet-module doesn't exist"
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



run_nlojet_module () { 

    cd $BASEDIR/nlojet-module

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
	./standalone output/weight_c.root
	./stand      output/weight_c.root
    else
	echo "nlojet grid code did not complete correctly" 
    fi
    
}



ARGS=$*



install_appl_grid   $ARGS
install_pdf         $ARGS
#install_mcfm        $ARGS
#run_mcfm        
install_nlojet         $ARGS
install_nlojet_module  $ARGS
run_nlojet_module 
