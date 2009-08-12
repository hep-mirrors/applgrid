AUTOMAKE_OPTIONS = foreign
SUBDIRS = bits

includedir= @includedir@/nlo++
include_HEADERS = dis2jet.h dis3jet.h dis4jet.h


INCLUDES = -I$(top_srcdir)/proc-dis -I$(top_srcdir)/nlo-amp -I$(top_srcdir)/nlo-core \
-I$(top_srcdir)/nlo-ps -I$(top_srcdir)/nlo-hep/hep-lorentz -I$(top_srcdir)/nlo-hep/hep-random \
-I$(top_srcdir)/nlo-hep/hep-container -I$(top_srcdir)/nlo-hep/hep-stat \
-I$(top_srcdir)/nlo-hep/hep-support

noinst_LTLIBRARIES = libdis-nlo++.la


libdis_nlo___la_SOURCES = event.cc phasespace.cc process.cc weight.cc disjet.cc \
dis2jet.cc dis3jet.cc dis4jet.cc
