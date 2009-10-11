##########################################################################
#
#   File:         makefile
#
#                 Copyright (C) M.Sutton (sutt@hep.ucl.ac.uk) 
#
#   Description:  build the appl_grid project
#                 
#
#   Created:    M.Sutton (sutt@hep.ucl.ac.uk) Fri Mar  7 08:29:08 CET 2008 
#   
#   $Id: makefile 20 2008-05-05 01:55:42Z sutt $                
#                   
######################################################################### 

ifdef BASEDIR
INSTALLBASE = $(BASEDIR)
else
INSTALLBASE = ../..
endif

ARCH_TYPE = $(shell setarch.sh)

AGBASE = ..
SRCDIR = .
INCDIR = ../appl_grid
# INCDIR = .
LIBDIR = ../lib/$(ARCH_TYPE)
OBJDIR = ../obj/$(ARCH_TYPE)
EXEDIR = ../exe/$(ARCH_TYPE)


# Root variables
ROOTCFLAGS   = $(shell root-config --cflags)
ROOTLIBS     = $(shell root-config --libs)
ROOTGLIBS    = $(shell root-config --glibs)

# CLHEPLIB    = $(shell clhep-config --libs)
# CLHEPCFLAGS = $(shell clhep-config --include)
CLHEPCFLAGS = 

# HOPPETFLAG = $(shell hoppet.sh $(BASEDIR)/include/hoppet_v1.h) -I$(BASEDIR)/include


# BASEDIR   = $(shell ./cwd.sh -h)
# BASEDIR   = .

# CXX = g++
# CC  = gcc
# F77 = g77

CSOURCES := $(wildcard *.cpp)
DEPENDENCIES = $(patsubst %.cpp, %.d, $(CSOURCES))


CFLAGS   += -O3 -fPIC -I$(AGBASE)
CXXFLAGS += -O3 -fPIC $(ROOTCFLAGS) $(CLHEPCFLAGS) -I$(AGBASE) $(HOPPETFLAG) $(HOPPETINCS) -Wno-deprecated  





GLIBS += $(ROOTGLIBS) $(ROOTLIBS)

LDFLAGS =

AR      = ar
ARFLAGS = -r

RANLIB   = ranlib
RANFLAGS = 

RM  = rm
TAR = tar


# CXXFLAGS += -pg $(shell root-config --cflags) -I. -I$(HOME)/include
CXXFLAGS += -O3 -pg $(shell root-config --cflags) -D_DIR_=\"$(PWD)\" 
# -I. 

F77 = g77
FFLAGS = 

AR = ar

RANLIB = ranlib

CP = cp

# FLIBS = -L/usr/local/lib -lg2c 
FLIBS = -L/sw/lib -lg2c 
CLIBS = -L/usr/lib/gcc/darwin/3.3 -lstdc++ 

RLIBS = $(shell root-config --libs)
GLIBS = $(shell root-config --glibs)

# LIBS = -lLHAPDF -L$(HOME)/lib -lreadcards -llcg -L. -lAPPLgrid $(FLIBS) $(RLIBS) 
LIBS = -lLHAPDF -L$(LIBDIR) -lAPPLgrid $(FLIBS) $(RLIBS) 

SRCFILES = $(FFILES)
OBJFILES = $(FFILES:.f=.o) main.o 
LOBJFILES = $(OBJDIR)/fappl_grid.o      \
	    $(OBJDIR)/appl_grid.o       \
	    $(OBJDIR)/appl_igrid.o      \
	    $(OBJDIR)/appl_pdf.o        \
	    $(OBJDIR)/nlojet_pdf.o      \
	    $(OBJDIR)/nlojetpp_pdf.o    \
	    $(OBJDIR)/jetrad_pdf.o      \
	    $(OBJDIR)/mcfmz_pdf.o       \
	    $(OBJDIR)/mcfmw_pdf.o       \
	    $(OBJDIR)/dis_pdf.o         \
	    $(OBJDIR)/TFileString.o     \
	    $(OBJDIR)/TFileStringDict.o \
	    $(OBJDIR)/appl_timer.o      \
            $(OBJDIR)/SparseMatrix3d.o  \
            $(OBJDIR)/fastnlo.o         \
	    $(OBJDIR)/fastnlo-dis.o     \
	    $(OBJDIR)/hoppet_init.o      



# HEADERS = APPLgrid.h APPLigrid.h
HEADERS = $(INCDIR)/appl_grid.h    \
  	  $(INCDIR)/appl_igrid.h   \
          $(INCDIR)/appl_pdf.h     \
          $(INCDIR)/appl_timer.h   \
	  $(INCDIR)/mcfmw_pdf.h    \
	  $(INCDIR)/mcfmz_pdf.h    \
	  $(INCDIR)/nlojet_pdf.h   \
	  $(INCDIR)/nlojetpp_pdf.h \
	  $(INCDIR)/jetrad_pdf.h   \
	  $(INCDIR)/dis_pdf.h      \
	  $(INCDIR)/SparseMatrix3d.h \
	  $(INCDIR)/axis.h           \
	  $(INCDIR)/sparse.h       \
	  $(INCDIR)/Directory.h    \
	  $(INCDIR)/tsparse_base.h \
	  $(INCDIR)/tsparse1d.h    \
	  $(INCDIR)/tsparse2d.h    \
	  $(INCDIR)/tsparse3d.h    \
	  $(INCDIR)/fastnlo.h      \
	  $(INCDIR)/fastnlo-dis.h  \
	  $(INCDIR)/hoppet_init.h 


.SUFFIXES : .c .f .cxx 

# compile
$(OBJDIR)/%.o : $(SRCDIR)/%.cxx $(HEADERS)
	$(CXX) $(CXXFLAGS) -c $<  -o $@

$(OBJDIR)/%.o : $(SRCDIR)/%.c $(HEADERS)
	$(CC) $(CFLAGS) -c $<  -o $@


$(OBJECTS) : $(HEADERS)

.PHONY : 



install : mkdirs installheaders installlibs

installheaders : $(HEADERS)
	./mkdir.sh  $(INSTALLBASE)/include
	./mkdir.sh  $(INSTALLBASE)/include/appl_grid
	cp $(HEADERS) $(INSTALLBASE)/include/appl_grid    

installlibs : lib timerlib
	./mkdir.sh  $(INSTALLBASE)/lib
	cp $(LIBDIR)/lib* $(INSTALLBASE)/lib

lib : $(LIBDIR)/libAPPLgrid.a 

# $(LIBDIR)/libgtttimer.a

$(LIBDIR)/libAPPLgrid.a : $(LOBJFILES)
	$(AR) -r $@ $(LOBJFILES)
	$(RANLIB) $@


# main.o : main.cxx
#	$(CXX) -c main.cxx

$(NAME): $(OBJFILES)
#	 $(FC) -o $@ $(FFLAGS) $(OBJFILES) -lLHAPDF -L$(HOME)/lib -lreadcards
	 $(CXX) -o $@ $(FFLAGS) $(OBJFILES) $(LIBS)

# tester : main.o install
#	$(CXX) -pg -o $@ main.o $(CFFLAGS) $(LIBS)

tester : $(OBJDIR)/test.o
	$(CXX) -pg -o $@ $(OBJDIR)/test.o $(CFFLAGS) $(LIBS)

crss : $(OBJDIR)/crss.o $(OBJDIR)/genpdf.o install
	$(CXX) -pg -o $@ $(OBJDIR)/crss.o $(OBJDIR)/genpdf.o $(CFFLAGS) $(LIBS)


toaster : tester
	$(CP) $< $@


# grid : $(OBJDIR)/main.o $(OBJDIR)/fmain.o lib
#	$(CXX) -o $@ $(OBJDIR)/main.o $(OBJDIR)/fmain.o $(LIBS)

grid : $(OBJDIR)/main.o lib
	$(CXX) -o $@ $(OBJDIR)/main.o  $(LIBS) -L$(HOME)/lib -ldglap -lgfortran



timerlib : 

# $(LIBDIR)/libgtttimer.a

# $(LIBDIR)/libgtttimer.a : $(OBJDIR)/gtttimer.o
# 	ar -r $@ $<
#	ranlib $@


TFileStringDict.cxx : $(INCDIR)/TFileString.h TFileString.cxx 
	rootcint -f $@ -c $(INCDIR)/TFileString.h 


timerinstall : 

# $(LIBDIR)/libgtttimer.a 


dirs : mkdirs
mkdirs : 
	./createdirs.sh


clean: 
	$(RM) -f $(SRCDIR)/*Dict.{h,cxx}
	$(RM) -f $(OBJDIR)/*.o $(LIBDIR)/*.so $(LIBDIR)/lib*.a 


tidy:
	$(RM) -f *~

distclean: clean tidy
	$(RM) -rf ../{obj,lib,exe}

# 	$(RM) $(SRCDIR)/appl_grid	


archive: clean tidy
	cd ../.. ; tar -czf appl_grid.tgz appl_grid/src


