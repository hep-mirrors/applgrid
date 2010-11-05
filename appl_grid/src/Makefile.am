##########################################################################
#
#   File:         Makefile.am
#
#                 Copyright (C) M.Sutton (sutt@cern.ch) 
#
#   Description:  generate the Makefile for the appl_grid project
#                 
#
#   Created:    M.Sutton (sutt@hep.ucl.ac.uk) Fri  5 Nov 2010 00:05:29 GMT
#   
#   $Id: makefile 20 2008-05-05 01:55:42Z sutt $                
#                   
#########################################################################

 
lib_LTLIBRARIES= libAPPLgrid.la libfAPPLgrid.la 
#bin_PROGRAMS = 
libAPPLgrid_la_SOURCES = \
	appl_grid.cxx		appl_igrid.cxx  	appl_pdf.cxx 		appl_timer.cxx \
	fastnlo.cxx 		SparseMatrix3d.cxx	hoppet_init.cxx	\
	nlojet_pdf.cxx		nlojetpp_pdf.cxx 	mcfmz_pdf.cxx 		mcfmw_pdf.cxx  \
	jetrad_pdf.cxx		dis_pdf.cxx		TFileString.cxx		TFileStringDict.cxx

#	fappl_grid.cxx

#bin_PROGRAMS = 
libfAPPLgrid_la_SOURCES = fappl_grid.cxx


library_includedir=$(includedir)/appl_grid
library_include_HEADERS=../appl_grid/*.h
ROOTCFLAGS= $(shell root-config --cflags)
ROOTLIBS= $(shell root-config --nonew --libs)
ROOTGLIBS= $(shell root-config --nonew --glibs)
LHAPDFDIR= $(shell lhapdf-config --prefix)
HOPPETDIR= $(shell hoppet-config --prefix)
HOPPETCFLAGS= $(shell hoppet-config --cxxflags)

ROOTARCH = $(findstring -m64, $(ROOTCFLAGS) )

# GSLLIBS=$(shell gsl-config --libs)
# WDIR = $(shell pwd)
ROOTCINT=rootcint
AM_CXXFLAGS= $(ROOTARCH) -O2 -Wall -fPIC  -I../appl_grid $(ROOTCFLAGS) -I$(LHAPDFDIR)/include $(HOPPETCFLAGS) 
AM_CFLAGS=   $(ROOTARCH) -O2 -Wall -fPIC  -I../appl_grid $(ROOTCFLAGS) -I$(LHAPDFDIR)/include $(HOPPETCFLAGS) 
#INCLUDES=-I$(LHAPDFDIR)/include
# AM_LDFLAGS= -O -L$(LHAPDFDIR)/lib -lLHAPDF  $(ROOTLIBS) $(GSLLIBS) $(FCLIBS)
AM_LDFLAGS= $(ROOTARCH) -O -L$(LHAPDFDIR)/lib -lLHAPDF  $(ROOTLIBS) $(FCLIBS)
# program_LDFLAGS= $(AM_LDFLAGS) 
AM_FCFLAGS= -c
AM_FCLIBS= -lg2c
AM_SOFLAGS= -shared
CINT= rootcint
clean-local:
	rm *.o *.lo *Dict*
%.o : %.cxx 
	@ echo ROOTARCH $(ROOTARCH) 
	$(CXX) $(AM_CXXFLAGS) -c $<
%.o : %.c 
	$(CC) $(AM_CFLAGS) -c $<
%Dict.cxx : ../appl_grid/%.h %.cxx
	$(CINT) -f $@ -c $< 

#../appl_grid/$*LinkDef.h
