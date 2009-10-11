AUTOMAKE_OPTIONS = foreign
SUBDIRS = bits

includedir= @includedir@/nlo++
include_HEADERS = photo2jet.h photo3jet.h photo4jet.h


INCLUDES = -I$(top_srcdir)/proc-photo -I$(top_srcdir)/proc-hhc \
-I$(top_srcdir)/nlo-amp -I$(top_srcdir)/nlo-ps -I$(top_srcdir)/nlo-core \
-I$(top_srcdir)/nlo-hep/hep-lorentz -I$(top_srcdir)/nlo-hep/hep-random \
-I$(top_srcdir)/nlo-hep/hep-container -I$(top_srcdir)/nlo-hep/hep-stat \
-I$(top_srcdir)/nlo-hep/hep-support


noinst_LTLIBRARIES = libphoto-nlo++.la


libphoto_nlo___la_SOURCES = phasespace.cc process.cc weight.cc \
photojet.cc photo2jet.cc photo3jet.cc photo4jet.cc

