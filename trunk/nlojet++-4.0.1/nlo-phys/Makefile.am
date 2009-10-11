AUTOMAKE_OPTIONS = foreign
SUBDIRS = bits

includedir= @includedir@/nlo++

include_HEADERS = kT_clus.h cone_clus.h


INCLUDES = -I$(top_srcdir)/nlo-phys -I$(top_srcdir)/nlo-core \
-I$(top_srcdir)/nlo-hep/hep-lorentz -I$(top_srcdir)/nlo-hep/hep-random \
-I$(top_srcdir)/nlo-hep/hep-container -I$(top_srcdir)/nlo-hep/hep-stat \
-I$(top_srcdir)/nlo-hep/hep-support -I$(top_srcdir)/proc-epa -I$(top_srcdir)/proc-dis \
-I$(top_srcdir)/proc-hhc -I$(top_srcdir)/proc-hhc2ph 

noinst_LTLIBRARIES = libphys-nlo++.la

libphys_nlo___la_SOURCES = cone_seedless.cc kT_clus_epa.cc \
kT_clus_ini.cc kT_clus_dis.cc kT_clus_long.cc 

