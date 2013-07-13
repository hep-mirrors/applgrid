EXTRA_DIST = applgrid-config.in 
bin_SCRIPTS = applgrid-config


bin_PROGRAMS = applgrid-combine
applgrid_combine_SOURCES = combine.cxx

ROOTARCH=$(shell root-config --ldflags)
ROOTCFLAGS= $(shell root-config --cflags)
ROOTLIBS= $(shell root-config --nonew --libs)
ROOTGLIBS= $(shell root-config --nonew --glibs)

APPLCFLAGS= -I../src -I../appl_grid
APPLLIBS= -L../src -lAPPLgrid

AM_LDFLAGS= $(ROOTARCH) -O2  $(ROOTLIBS) $(APPLLIBS)
AM_CXXFLAGS= -O2 -Wall -fPIC  -I. -I.. -I$(srcdir)  $(ROOTCFLAGS) $(APPLCFLAGS)


clean-local:
	rm -rf *.o *.lo 
.cxx.o :   
	$(CXX) $(AM_CXXFLAGS) -c $<

