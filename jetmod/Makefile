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
# CXXFLAGS += -I$(INSTALLBASE)/include/nlo++ -I$(INSTALLBASE)/include
CXXFLAGS += -I$(INSTALLBASE)/include/nlo++ -I$(INSTALLBASE)/include


HOPPETINCS =  $(shell hoppet-config --cxxflags)
HOPPETLIBS =  $(shell hoppet-config --libs)

CXXFLAGS += $(HOPPETINCS)



# ROOT	
# ROOTINCS = -I$(shell root-config --incdir)

ROOTINCS = -I$(shell root-config --incdir)
ROOTLIBS = $(shell root-config --libs)

# F2CFLAG = -lg2c
#  -L/usr/lib/gcc/i386-redhat-linux5E/4.1.2/ -lgfortran
ROOTARCH = $(findstring -m64, $(ROOTINCS) )

FFLAGS += -g $(ROOTARCH)


CXXFLAGS += $(ROOTARCH) $(ROOTINCS)

#LHAPDF
LHAPDFDIR = $(shell lhapdf-config --prefix)
LHAPDFINCS = -I$(LHAPDFDIR)/include
LHAPDFLIBS = -L$(LHAPDFDIR)/lib -lLHAPDF

APPLCXXFLAGS = $(shell applgrid-config --cxxflags)
APPLLIBS     = $(shell applgrid-config --ldcflags)

CXXFLAGS += $(APPLCXXFLAGS) $(LHAPDFINCS)

LD_FLAGS += $(ROOTARCH) -L$(LHAPDFDIR)/lib -lLHAPDF
LD_FLAGS += $(APPLLIBS)

LDFLAGS += $(ROOTARCH) 
#
#F90 libraries
#
# F90LIB= -L`which lf95 | sed s:bin/lf95:lib:` -lfj9i6 -lfj9f6 -lfj9e6 -lfccx86_6a
FRTLIB   = $(shell gfortran $(CXXFLAGS) -print-file-name=libgfortran.a)
FRTLLIB  = -L$(subst /libgfortran.a, ,$(FRTLIB) ) -lgfortran 
# FRTLLIB = 
# FRTLLIB = $(shell gfortran.sh )
# FRTLLIB = -L/usr/local/gfortran/lib/x86_64  -lgfortran
# FRTLLIB = -lg2c -lgfortran
# FRTLLIB = -lgfortran
F90LIB = $(FRTLLIB)
#PDFCONVLIBS = -L$(INSTPATH)/pdf-conv-1.0/src -ldglap $(F90LIB)
#PDFCONVLIBS = -L$(INSTALLBASE)/lib -ldglap $(F90LIB) $(LHAPDFLIBS)

PDFCONVLIBS = $(HOPPETLIBS) $(F90LIB) $(LHAPDFLIBS)


# macros to setup contact interactions and reference histograms
LHA = -DUSELHAPDF
REF = -DREN_REFERENCE

ALLMACRO  = $(LHA) $(REF)
ALLMACRO1 = $(LHA),$(REF)

ALL_LIBS = -L$(INSTALLBASE)/lib,-L$(INSTALLBASE)/lib,-lAPPLgrid,-L$(LHAPDFDIR)/lib,-lLHAPDF,`root-config --libs | sed 's/ /,/g'`



NLOCXXFLAGS = -Wx,$(ARCH),$(ALLMACRO1),-O3,-Wall,$(ROOTINCS),-I$(INSTALLBASE)/include,-I$(INSTALLBASE)/include/lhpdf,$(APPLCXXFLAGS)
NLOFFLAGS   = -Wf,$(ARCH),-O3,-Wall
NLOLIBFLAGS = -Wl,$(ARCH),-lpthread,$(ALL_LIBS),-L$(INSTALLBASE)/lib,-lAPPLgrid,nlogrid.lo


NLOINC = -I$(INSTALLBASE)/include/nlo++ 

install : all
all : fillgrid stand


standalone: standalone.o  InitMyPdf.o mypdfevolve.o
	$(CXX) -o standalone $(LDFLAGS) $(ALLMACRO)  standalone.o InitMyPdf.o mypdfevolve.o    $(APPLLIBS) $(ROOTLIBS) $(PDFCONVLIBS) 

stand: stand.o  InitMyPdf.o mypdfevolve.o
	@ echo FRTLIB $(FRTLIB) $(FRTLLIB)
	@ echo ROOTARCH $(ROOTARCH)
	$(CXX) -o stand $(LDFLAGS) $(ALLMACRO)  stand.o InitMyPdf.o mypdfevolve.o    $(APPLLIBS) $(ROOTLIBS) $(PDFCONVLIBS) 


fillgrid: LHAPDF.h nlogrid.o nlogrid.lo
	$(INSTALLBASE)/bin/create-nlojet-user lhapdf_string_interface.f fillgrid.cxx -o fillgrid $(NLOCXXFLAGS) $(NLOFFLAGS) $(NLOLIBFLAGS) 


.SUFFIXES : .cxx .o .f .c

.f.o :
	$(F77) $(FFLAGS) -c $<

.cxx.o:	 
	$(CXX) -c $< $(CXXFLAGS) $(ALLMACRO) -I$(INSTALLBASE)/include -I$(INSTALLBASE)/include/lhpdf $(APPLCXXFLAGS)



# obsolete single compilation rules 

# nlogrid.o: nlogrid.cxx
# 	$(CXX) $(ALLMACRO) -g -c $(CXXFLAGS) nlogrid.cxx $(ROOTINCS) $(NLOINC) $(APPLCXXFLAGS) 

# standalone.o:	standalone.cxx 
#	$(CXX) -c $(CXXFLAGS) $(ALLMACRO) standalone.cxx -I$(INSTALLBASE)/include -I$(INSTALLBASE)/include/lhpdf $(APPLCXXFLAGS)

NLOLIBTOOLFLAGS = -Wc,$(ARCH),$(ROOTINCS),-I$(INSTALLBASE)/include,-I$(INSTALLBASE)/include/nlo++,-I$(INSTALLBASE)/include/lhpdf,$(APPLCXXFLAGS)

nlogrid.lo: nlogrid.cxx 
#	$(INSTPATH)/libexec/nlo++/libtool-nlojet++ --mode=compile $(CXX) $(ALLMACRO) -c nlogrid.cxx $(NLOLIBTOOLFLAGS)
	$(INSTALLBASE)/libexec/nlo++/libtool-nlojet++ --mode=compile $(CXX) $(ALLMACRO) $(CXXFLAGS) -c nlogrid.cxx $(NLOLIBTOOLFLAGS)



clean:
	rm -rf ./.libs ./.obj *.lo *.o *.la standalone
