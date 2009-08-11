AUTOMAKE_OPTIONS = foreign
SUBDIRS = bits

includedir= @includedir@/nlo++
noinst_HEADERS = defmacros.h
include_HEADERS = ampq2g1l2.h ampq2g2l2.h ampq2g3l2.h ampq4l2.h ampq4g1l2.h \
ampg4.h ampq2g2.h ampq4.h \
ampg5.h ampq2g3.h ampq4g1.h \
ampg6.h ampq2g4.h ampq4g2.h ampq6.h \
ampq2g1p2.h ampq2g2p2.h ampq4p2.h \
ampq2g1p1.h ampq2g2p1.h ampq4p1.h ampq2g3p1.h ampq4g1p1.h \
ampp6.h


INCLUDES = -I$(top_srcdir)/nlo-amp -I$(top_srcdir)/nlo-core \
-I$(top_srcdir)/nlo-hep/hep-lorentz -I$(top_srcdir)/nlo-hep/hep-random \
-I$(top_srcdir)/nlo-hep/hep-container -I$(top_srcdir)/nlo-hep/hep-stat \
-I$(top_srcdir)/nlo-hep/hep-support

noinst_LTLIBRARIES = libamp-nlo++.la


libamp_nlo___la_SOURCES = li2.cc cl2.cc ampbase.cc \
ampq2g1l2.cc ampq2g2l2.cc ampq2g3l2.cc ampq4l2.cc ampq4g1l2.cc \
ampg4.cc ampq2g2.cc ampq4.cc \
ampg5.cc ampq2g3.cc ampq4g1.cc \
ampg6.cc  ampq2g4.cc ampq4g2.cc ampq6.cc \
ampq2g1p2.cc ampq2g2p2.cc ampq4p2.cc \
ampq2g1p1.cc ampq2g2p1.cc ampq4p1.cc ampq2g3p1.cc ampq4g1p1.cc \
ampp6.cc
