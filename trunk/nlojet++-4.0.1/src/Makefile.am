AUTOMAKE_OPTIONS = foreign

includedir= @includedir@/nlo++
include_HEADERS = nlo++-module_add.h

INCLUDES = -I$(top_srcdir)/nlo-amp -I$(top_srcdir)/nlo-ps \
-I$(top_srcdir)/nlo-core -I$(top_srcdir)/nlo-hep/hep-lorentz -I$(top_srcdir)/nlo-hep/hep-random \
-I$(top_srcdir)/nlo-hep/hep-container -I$(top_srcdir)/nlo-hep/hep-stat \
-I$(top_srcdir)/nlo-hep/hep-support  \
-I$(top_srcdir)/proc-epa -I$(top_srcdir)/proc-dis \
-I$(top_srcdir)/proc-hhc -I$(top_srcdir)/proc-hhc2ph \
-I$(top_srcdir)/proc-photo -I$(top_srcdir)/libltdl -I../nlo-core


noinst_HEADERS = nlojet++.h

bin_PROGRAMS = nlojet++
bin_SCRIPTS = create-nlojet-user

processes = main_epa.cc main_dis.cc main_hhc.cc main_hhc2ph.cc main_photo.cc
nlojet___SOURCES = $(processes) nlojet++.cc main_calc.cc main_add.cc main.cc 
nlojet___LDADD = ../libnlo++.la $(LIBLTDL)

