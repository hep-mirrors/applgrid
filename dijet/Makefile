#Set the Install path as the current directory
# need more standardised directory structure and setup
#INSTPATH=`pwd`
# GPS alternative for when one has a standard nlojet (v4) installation

# BASE   = $(HOME)

# INSTALLBASE = $(BASEDIR)

NLODIR = $(INSTALLBASE)

INSTPATH= `which nlojet++ | sed s:/bin/nlojet\+\+::`
# CXX = g++
F77 = gfortran
CXXFLAGS += -O3 -fPIC -g 
CXXFLAGS += -I$(INSTALLBASE)/include/nlo++ -I$(INSTALLBASE)/include


HOPPETINCS =  $(shell hoppet-config --cxxflags)
HOPPETLIBS =  $(shell hoppet-config --libs)

CXXFLAGS += $(HOPPETINCS)


# ROOT	

ROOTINCS = -I$(shell root-config --incdir)
ROOTLIBS = $(shell root-config --libs)
ROOTCFLAGS = $(shell root-config --cflags)
ROOTARCH = $(findstring -m64, $(ROOTCFLAGS) )

FFLAGS += -g $(ROOTARCH)


CXXFLAGS += $(ROOTARCH) $(ROOTINCS)

# LHAPDF
LHAPDFDIR = $(shell lhapdf-config --prefix)
LHAPDFINCS = -I$(LHAPDFDIR)/include
LHAPDFLIBS = -L$(LHAPDFDIR)/lib -lLHAPDF

# applgrid
APPLCXXFLAGS = $(shell applgrid-config --cxxflags)
APPLLIBS     = $(shell applgrid-config --ldcflags)

# fastjet
FASTJETINCS = $(shell fastjet-config --cxxflags)
FASTJETLIBS = $(shell fastjet-config --libs)
FASTJETFORNLO = -L$(shell fastjet-config --prefix)/lib,-lfastjet


CXXFLAGS += $(APPLCXXFLAGS) $(LHAPDFINCS) $(FASTJETINCS)

LD_FLAGS += $(ROOTARCH) -L$(LHAPDFDIR)/lib -lLHAPDF
LD_FLAGS += $(APPLLIBS)

LDFLAGS += $(ROOTARCH) 


#F90 libraries
# F2CFLAG = -lg2c
#  -L/usr/lib/gcc/i386-redhat-linux5E/4.1.2/ -lgfortran
# F90LIB= -L`which lf95 | sed s:bin/lf95:lib:` -lfj9i6 -lfj9f6 -lfj9e6 -lfccx86_6a
FRTLIB   = $(shell gfortran $(CXXFLAGS) -print-file-name=libgfortran.a)
FRTLLIB  = -L$(subst /libgfortran.a, ,$(FRTLIB) ) -lgfortran 
# FRTLLIB = $(shell gfortran.sh )
F90LIB = $(FRTLLIB)
#PDFCONVLIBS = -L$(INSTPATH)/pdf-conv-1.0/src -ldglap $(F90LIB)
#PDFCONVLIBS = -L$(INSTALLBASE)/lib -ldglap $(F90LIB) $(LHAPDFLIBS)

PDFCONVLIBS = $(HOPPETLIBS) $(F90LIB) $(LHAPDFLIBS)


# macros to setup contact interactions and reference histograms
LHA = -DUSELHAPDF
REF = -DREN_REFERENCE

ALLMACRO  = $(LHA) $(REF)
ALLMACRO1 = $(LHA),$(REF)

ALL_LIBS = -L$(INSTALLBASE)/lib,$(FASTJETFORNLO),-lAPPLgrid,-L$(LHAPDFDIR)/lib,-lLHAPDF,`root-config --libs | sed 's/ /,/g'`



NLOCXXFLAGS = -Wx,$(ROOTARCH),$(ALLMACRO1),-O3,-Wall,$(ROOTARCH),$(ROOTINCS),$(FASTJETINCS),-I$(INSTALLBASE)/include,-I$(INSTALLBASE)/include/lhpdf,$(APPLCXXFLAGS)
NLOFFLAGS   = -Wf,$(ROOTARCH),-O3,-Wall,$(FASTJETINCS)
NLOLIBFLAGS = -Wl,$(ROOTARCH),-lpthread,$(ALL_LIBS),-L$(INSTALLBASE)/lib,-lAPPLgrid,nlogrid.lo,fjClustering.lo

NLOINC = -I$(INSTALLBASE)/include/nlo++ 

NLOLIBTOOLFLAGS = -Wc,$(ROOTARCH),$(ROOTINCS),-I$(INSTALLBASE)/include,-I$(INSTALLBASE)/include/nlo++,-I$(INSTALLBASE)/include/lhpdf,$(APPLCXXFLAGS)



install : all
all : fillgrid standalone



standalone: standalone.o  InitMyPdf.o mypdfevolve.o
	$(CXX) -o standalone $(LDFLAGS) $(ALLMACRO) $(HOPPETINCS) standalone.o InitMyPdf.o mypdfevolve.o    $(APPLLIBS) $(ROOTLIBS) $(PDFCONVLIBS) $(HOPPETLIB)

standaloneScalePDF: standaloneScalePDF.o  InitMyPdf.o mypdfevolve.o
	$(CXX) -o standaloneScalePDF $(LDFLAGS) $(ALLMACRO) $(HOPPETINCS) standaloneScalePDF.o InitMyPdf.o mypdfevolve.o    $(APPLLIBS) $(ROOTLIBS) $(PDFCONVLIBS) $(HOPPETLIB)

fillgrid: LHAPDF.h nlogrid.o nlogrid.lo fjClustering.o fjClustering.lo
	$(INSTALLBASE)/bin/create-nlojet-user lhapdf_string_interface.f fillgrid.cpp -o fillgrid $(NLOCXXFLAGS) $(NLOFFLAGS) $(NLOLIBFLAGS) 


.SUFFIXES : .cxx .o .f .c .lo

.f.o :
	$(F77) $(FFLAGS) -c $<

.cpp.o:	 
	$(CXX) -c $< $(CXXFLAGS) $(ALLMACRO) -I$(INSTALLBASE)/include -I$(INSTALLBASE)/include/lhpdf $(APPLCXXFLAGS)

.cpp.lo:  
	$(INSTALLBASE)/libexec/nlo++/libtool-nlojet++ --mode=compile $(CXX) $(ALLMACRO) $(CXXFLAGS) -c $< $(NLOLIBTOOLFLAGS)


clean:
	rm -rf ./.libs ./.obj *.lo *.o *.la standalone
