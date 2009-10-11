AUTOMAKE_OPTIONS = foreign
SUBDIRS = bits

includedir= @includedir@/nlo++

include_HEADERS = distribution.h statistic.h 


INCLUDES = -I$(top_srcdir)/nlo-hep/hep-stat -I$(top_srcdir)/nlo-hep/hep-support \
-I$(top_srcdir)/nlo-core
noinst_LTLIBRARIES = libhep-stat.la

libhep_stat_la_SOURCES = sample_traits.cc statistic.cc


