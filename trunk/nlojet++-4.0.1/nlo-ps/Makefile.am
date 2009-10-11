AUTOMAKE_OPTIONS = foreign
SUBDIRS = bits

includedir= @includedir@/nlo++

include_HEADERS =

INCLUDES = -I$(top_srcdir)/nlo-ps -I$(top_srcdir)/nlo-core \
-I$(top_srcdir)/nlo-hep/hep-lorentz -I$(top_srcdir)/nlo-hep/hep-random \
-I$(top_srcdir)/nlo-hep/hep-container -I$(top_srcdir)/nlo-hep/hep-stat \
-I$(top_srcdir)/nlo-hep/hep-support

noinst_LTLIBRARIES = libps-nlo++.la

libps_nlo___la_SOURCES = rambo.cc dipole_emission.cc \
phasespace_n0i0f0.cc phasespace_n1i1f0.cc phasespace_n0i2f0.cc \
phasespace_n2i2f0.cc
