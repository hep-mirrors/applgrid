AC_INIT(APPLgrid, 1.0.3, sutt@cern.ch)
AM_INIT_AUTOMAKE
AM_CONFIG_HEADER(amconfig.h:amconfig.in)
#AC_PROG_FC
#AC_CHECK_PROGS(FC, [gfortran g77 f77 fort77 f90], no)
#AC_PROG_FC
#AC_LANG_PUSH(Fortran)
#AC_FC_SRCEXT(f90)
#AC_FC_LIBRARY_LDFLAGS
#AC_FC_WRAPPERS
#AC_FC_FREEFORM
#AC_FC_DUMMY_MAIN
#AC_LANG_POP(Fortran)
AC_PROG_CC
AC_PROG_CXX 
AM_PROG_LIBTOOL
AC_HEADER_STDC
AC_FUNC_UTIME_NULL
## Set Fortran compiler behaviour
#if test "x$FCFLAGS" == "x"; then 
#  FCFLAGS="-O2"
#fi
## Try to respect users' Fortran compiler variables
#if test "x$FC" == "x"; then 
#  if test "x$F77" == "x"; then 
#FC="$GFORTRAN"
#FC=gfortran
#  else
#    FC="$F77"
#  fi
#fi
AC_PATH_PROG(ROOTPATH,root,no)
if test "$ROOTPATH" = "no"; then
    AC_MSG_WARN([root was not found - needed for applgrid storage])
fi

# lhapdf isn't technically needed but leav it in anyhow
AC_PATH_PROG(LHAPDFPATH,lhapdf-config,no)

AC_PATH_PROG(HOPPETPATH,hoppet-config,no)
if test "$HOPPETPATH" = "no"; then
    AC_MSG_WARN([HOPPET was not found - needed for arbitrary applgrid scale evolution])
else 
    CPPFLAGS+=" -DHOPPET "
    CXXFLAGS+=" -DHOPPET "
fi


AC_CONFIG_FILES([bin/Makefile bin/applgrid-config])


AC_OUTPUT(Makefile src/Makefile)
