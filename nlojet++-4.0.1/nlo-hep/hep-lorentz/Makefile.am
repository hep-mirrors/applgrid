AUTOMAKE_OPTIONS = foreign
SUBDIRS = bits

includedir= @includedir@/nlo++

include_HEADERS = lorentzvector.h threevector.h 


INCLUDES = -I$(top_srcdir)/nlo-hep/hep-lorentz
noinst_LTLIBRARIES = libhep-lorentz.la


libhep_lorentz_la_SOURCES = threevector.cc lorentzvector.cc

