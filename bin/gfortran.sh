#!/bin/csh -f

if ( $?ARCH ) then   
  set gf = `gfortran $ARCH -print-file-name=libgfortran.a` 
  echo "-L${gf:h} -lgfortran" 
endif


