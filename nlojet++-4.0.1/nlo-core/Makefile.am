AUTOMAKE_OPTIONS = foreign
SUBDIRS = bits

includedir= @includedir@/nlo++

include_HEADERS = 


INCLUDES = -I$(top_srcdir)/nlo-core -I$(top_srcdir)/nlo-hep/hep-lorentz -I$(top_srcdir)/nlo-hep/hep-random \
-I$(top_srcdir)/nlo-hep/hep-container -I$(top_srcdir)/nlo-hep/hep-stat -I$(top_srcdir)/nlo-hep/hep-support
noinst_LTLIBRARIES = libcore-nlo++.la


libcore_nlo___la_SOURCES = split.cc color.cc flkern.cc weight.cc \
basic_user.cc basic_user_result.cc
