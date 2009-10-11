AUTOMAKE_OPTIONS = foreign
SUBDIRS = bits

includedir= @includedir@/nlo++
include_HEADERS = hhc2ph1jet.h  hhc2ph2jet.h


INCLUDES = -I$(top_srcdir)/proc-hhc2ph -I$(top_srcdir)/nlo-amp \
-I$(top_srcdir)/nlo-core -I$(top_srcdir)/nlo-hep/hep-lorentz -I$(top_srcdir)/nlo-hep/hep-random \
-I$(top_srcdir)/nlo-hep/hep-container -I$(top_srcdir)/nlo-hep/hep-stat \
-I$(top_srcdir)/nlo-hep/hep-support

noinst_LTLIBRARIES = libhhc2ph-nlo++.la


libhhc2ph_nlo___la_SOURCES = process.cc weight.cc hhc2phjet.cc \
hhc2ph1jet.cc hhc2ph2jet.cc
