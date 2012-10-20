#Set the Install path as the current directory
# need more standardised directory structure and setup
#INSTPATH=`pwd`
# GPS alternative for when one has a standard nlojet (v4) installation

# BASE   = $(HOME)
# INSTALLBASE = $(BASEDIR)
NLODIR = $(INSTALLBASE)
INSTPATH= `which nlojet++ | sed s:/bin/nlojet\+\+::`

F77 = gfortran
FFLAGS     += -g 
#
#F90 libraries
#
FRTLLIB     = $(shell gfortran.sh )
F90LIB      = $(FRTLLIB)
#PDFCONVLIBS = $(HOPPETLIBS) $(F90LIB) $(LHAPDFLIBS)
PDFCONVLIBS = $(F90LIB) 

CXXFLAGS += -O3 -fPIC -g 
CXXFLAGS += -I$(INSTALLBASE)/include/nlo++ -I$(INSTALLBASE)/include


HOPPETINCS =  $(shell hoppet-config --cxxflags)
HOPPETLIBS =  $(shell hoppet-config --libs)
CXXFLAGS += $(HOPPETINCS)

ROOTINCS = -I$(shell root-config --incdir)
ROOTLIBS = $(shell root-config --libs)
CXXFLAGS   += $(ROOTINCS)

#LHAPDF
LHAPDFINCS = -I$(shell lhapdf-config --incdir) 
LHAPDFDIR   = $(shell lhapdf-config --prefix)/lib
LHAPDFLIBS  = -L$(LHAPDFDIR) -lLHAPDF

APPLCXXFLAGS = $(shell applgrid-config --cxxflags)
APPLLIBS     = $(shell applgrid-config --ldcflags)
#APPLCXXFLAGS= -I$(INSTALLBASE)/include
#APPLLIBS    = -L$(INSTALLBASE)/lib -lAPPLgrid

CXXFLAGS += $(APPLCXXFLAGS) $(LHAPDFINCS)
LD_FLAGS +=  $(LHAPDFLIBS)
LD_FLAGS += $(APPLLIBS)

FASTJET     = $(shell fastjet-config --prefix)
FASTJETLIBS = -L$(FASTJET)/lib,-lfastjet
FASTJETINCS = $(shell fastjet-config --cxxflags)
CXXFLAGS   += $(FASTJETINCS)

#
# macros to setup contact interactions and reference histograms
LHA = -DUSELHAPDF
REF = -DREN_REFERENCE

ALLMACRO  = $(LHA) $(REF)
ALLMACRO1 = $(LHA),$(REF)

ALL_LIBS           = -L$(INSTALLBASE)/lib,-lAPPLgrid,-L$(LHAPDFDIR),-lLHAPDF,`root-config --libs | sed 's/ /,/g'`

#NLOLIBTOOLFLAGS    = -Wc,$(ARCH),$(ROOTINCS),-I$(INSTALLBASE)/include,-I$(INSTALLBASE)/include/nlo++,-I$(INSTALLBASE)/include/lhpdf,$(APPLCXXFLAGS)
#NLOCXXFLAGS        = -Wx,$(ARCH),$(ALLMACRO1),-O3,-Wall,$(ROOTINCS),-I$(INSTALLBASE)/include,-I$(INSTALLBASE)/include/lhpdf,$(APPLCXXFLAGS),$(FASTJETINCS) 
#NLOFFLAGS          = -Wf,$(ARCH),-O3,-Wall
#NLOLIBFLAGS        = -Wl,$(ARCH),-lpthread,$(ALL_LIBS),$(FASTJETLIBS),-L$(INSTALLBASE)/lib,-lAPPLgrid,nlogrid.lo
#NLOINC             = -I$(INSTALLBASE)/include/nlo++ 


NLOLIBTOOLFLAGS    = -Wc,$(ARCH),$(ROOTINCS),-I$(INSTALLBASE)/include/nlo++,-I$(INSTALLBASE)/include,$(APPLCXXFLAGS)
NLOCXXFLAGS        = -Wx,$(ARCH),$(ALLMACRO1),-O3,-Wall,$(ROOTINCS),-I$(INSTALLBASE)/include,$(LHAPDFINCS),$(APPLCXXFLAGS),$(FASTJETINCS) 
NLOFFLAGS          = -Wf,$(ARCH),-O3,-Wall
#NLOLIBFLAGS        = -Wl,$(ARCH),-lpthread,-L$(INSTALLBASE)/lib,$(APPLLIBS),$(LHAPDFLIBS),$(ROOTLIBS),$(FASTJETLIBS),nlogrid.lo,
NLOLIBFLAGS        = -Wl,$(ARCH),-lpthread,-L$(INSTALLBASE)/lib,-lAPPLgrid,-L$(LHAPDFDIR),-lLHAPDF,`root-config --libs | sed 's/ /,/g'`,$(FASTJETLIBS),nlogrid.lo
#                                                                                                   ugly, ugly but need to remove commas..
install : all
all : fillgrid stand


standalone: standalone.o  InitMyPdf.o mypdfevolve.o
	$(CXX) -o standalone $(LDFLAGS) $(ALLMACRO)  standalone.o InitMyPdf.o mypdfevolve.o $(APPLLIBS) $(ROOTLIBS) $(LHAPDFLIBS)

stand: stand.o  InitMyPdf.o mypdfevolve.o
	@ echo FRTLIB $(FRTLIB) $(FRTLLIB)
	@ echo ARCH $(ARCH)
	$(CXX) -o stand $(LDFLAGS) $(ALLMACRO)  stand.o InitMyPdf.o mypdfevolve.o $(APPLLIBS) $(ROOTLIBS) $(LHAPDFLIBS)


fillgrid: LHAPDF.h nlogrid.lo
	$(INSTALLBASE)/bin/create-nlojet-user lhapdf_string_interface.f fillgrid.cxx -o fillgrid $(NLOCXXFLAGS) $(NLOFFLAGS) $(NLOLIBFLAGS) $(NLOLIBTOOLFLAGS),$(FASTJETINCS)
 


.SUFFIXES : .cxx .o .f .c

.f.o :
	$(F77) $(FFLAGS) -c $<

.cxx.o:	 
	$(CXX) -c $< $(CXXFLAGS) $(ALLMACRO) -I$(INSTALLBASE)/include -I$(INSTALLBASE)/include/lhpdf $(APPLCXXFLAGS)

nlogrid.lo: nlogrid.cxx
	$(INSTALLBASE)/libexec/nlo++/libtool-nlojet++ --mode=compile $(CXX) $(ALLMACRO) $(CXXFLAGS) -c nlogrid.cxx $(NLOLIBTOOLFLAGS)

clean:
	rm -rf ./.libs ./.obj *.lo *.o *.la standalone
