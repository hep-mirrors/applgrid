AUTOMAKE_OPTIONS = foreign
SUBDIRS = bits

includedir= @includedir@/nlo++
include_HEADERS = epa3jet.h epa4jet.h epa5jet.h


INCLUDES = -I$(top_srcdir)/proc-epa -I$(top_srcdir)/nlo-amp \
-I$(top_srcdir)/nlo-ps -I$(top_srcdir)/nlo-core \
-I$(top_srcdir)/nlo-hep/hep-lorentz -I$(top_srcdir)/nlo-hep/hep-random \
-I$(top_srcdir)/nlo-hep/hep-container -I$(top_srcdir)/nlo-hep/hep-stat \
-I$(top_srcdir)/nlo-hep/hep-support

noinst_LTLIBRARIES = libepa-nlo++.la


libepa_nlo___la_SOURCES = epajet.cc epa3jet.cc epa4jet.cc epa5jet.cc

