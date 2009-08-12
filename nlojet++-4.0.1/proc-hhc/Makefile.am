AUTOMAKE_OPTIONS = foreign
SUBDIRS = bits

includedir= @includedir@/nlo++
include_HEADERS = hhc2jet.h hhc3jet.h hhc4jet.h


INCLUDES = -I$(top_srcdir)/proc-hhc -I$(top_srcdir)/nlo-amp \
-I$(top_srcdir)/nlo-core -I$(top_srcdir)/nlo-hep/hep-lorentz -I$(top_srcdir)/nlo-hep/hep-random \
-I$(top_srcdir)/nlo-hep/hep-container -I$(top_srcdir)/nlo-hep/hep-stat \
-I$(top_srcdir)/nlo-hep/hep-support

noinst_LTLIBRARIES = libhhc-nlo++.la


libhhc_nlo___la_SOURCES = process.cc weight.cc hhcjet.cc hhc2jet.cc \
hhc3jet.cc hhc4jet.cc
