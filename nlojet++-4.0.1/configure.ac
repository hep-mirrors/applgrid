# Process this file with autoconf to produce a configure script.
AC_INIT(nlojet++, 4.0.1, Zoltan.Nagy@cern.ch)
AM_INIT_AUTOMAKE(nlojet++,4.0.1)

LIB_VERSION="4:0:1"
AC_SUBST(LIB_VERSION)

AC_CONFIG_SRCDIR([version.cc])
AM_CONFIG_HEADER(config.h)
AC_PREFIX_DEFAULT(${HOME})

AC_LIBLTDL_CONVENIENCE
AC_SUBST(INCLTDL)
AC_SUBST(LIBLTDL)
AC_LIBTOOL_DLOPEN

# Checks for programs.
AM_PROG_LIBTOOL
AC_PROG_CC
AC_PROG_CPP
AC_PROG_CXX
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET

# Checks for libraries.
AC_CHECK_LIB(m, main)
AC_F77_LIBRARY_LDFLAGS

# Configure libltdl
AC_CONFIG_SUBDIRS(libltdl)

# Checks for header files.
AC_CHECK_HEADERS([fenv.h limits.h stdlib.h string.h sys/wait.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_C_INLINE
AC_TYPE_PID_T
AC_TYPE_SIZE_T

# Checks for library functions.
AC_HEADER_STDC
AC_FUNC_MALLOC
AC_FUNC_STAT
AC_CHECK_FUNCS([lgamma gamma])

# Output
AC_CONFIG_FILES([Makefile src/Makefile
nlo-hep/Makefile
nlo-hep/hep-lorentz/Makefile nlo-hep/hep-lorentz/bits/Makefile
nlo-hep/hep-random/Makefile nlo-hep/hep-random/bits/Makefile
nlo-hep/hep-container/Makefile nlo-hep/hep-container/bits/Makefile
nlo-hep/hep-support/Makefile nlo-hep/hep-support/bits/Makefile
nlo-hep/hep-stat/Makefile nlo-hep/hep-stat/bits/Makefile
nlo-core/Makefile nlo-core/bits/Makefile
nlo-ps/Makefile nlo-ps/bits/Makefile
nlo-phys/Makefile nlo-phys/bits/Makefile
nlo-amp/Makefile nlo-amp/bits/Makefile
proc-epa/Makefile proc-epa/bits/Makefile
proc-dis/Makefile proc-dis/bits/Makefile
proc-hhc/Makefile proc-hhc/bits/Makefile
proc-hhc2ph/Makefile proc-hhc2ph/bits/Makefile
proc-photo/Makefile proc-photo/bits/Makefile
src/create-nlojet-user
])
AC_OUTPUT
