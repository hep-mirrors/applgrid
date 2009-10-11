# Process this file with autoconf to produce a configure script.
AC_INIT(lhpdf, 1.0.0, Zoltan.Nagy@durham.ac.uk)
AM_INIT_AUTOMAKE(lhpdf, 1.0.0, Zoltan.Nagy@durham.ac.uk)
AC_CONFIG_SRCDIR([lhpdf/lhpdf.h])
AM_CONFIG_HEADER([config.h])
AC_PREFIX_DEFAULT(${HOME})

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_CPP
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET
AC_PROG_RANLIB
AM_PROG_LIBTOOL
AC_PROG_F77

# Checks for libraries.
AC_F77_LIBRARY_LDFLAGS

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([stdlib.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST

# Checks for library functions.
AC_FUNC_MALLOC
AC_CHECK_FUNCS([pow sqrt])

AC_CONFIG_FILES([Makefile \
		 lhpdf/Makefile \
		 lhpdf++/Makefile \
		 test/Makefile \
])
AC_OUTPUT
