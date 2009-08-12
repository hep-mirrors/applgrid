AUTOMAKE_OPTIONS = foreign
SUBDIRS = bits

includedir= @includedir@/nlo++

include_HEADERS = random.h


INCLUDES = -I$(top_srcdir)/nlo-hep/hep-random
noinst_LTLIBRARIES = libhep-random.la


libhep_random_la_SOURCES = cmrg.cc mt.cc ranlux.cc ranlxd.cc ranlxs.cc gamma.cc \
gaussian.cc sphere.cc

