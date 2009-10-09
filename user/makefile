#Set the Install path as the current directory
# need more standardised directory structure and setup
# GPS alternative for when one has a standard nlojet (v4) installation


INSTALLBASE = $(BASEDIR)

CXX = g++
F77 = gfortran

FFLAGS += -g 

CXXFLAGS +=  -O3 -fPIC -g 
CXXFLAGS += -I$(INSTALLBASE)/include



ROOTINCS = -I$(shell root-config --incdir)
APPLCXXFLAGS = -I$(INSTALLBASE)/include

CXXFLAGS += $(ROOTINCS) $(APPLCXXFLAGS)



#LHAPDF
ROOTLIBS   = $(shell root-config --libs)
LHAPDFDIR  = $(shell lhapdf-config --prefix)/lib
LHAPDFLIBS = -L$(LHAPDFDIR) -lLHAPDF
APPLLIBS   = -L$(INSTALLBASE)/lib -lAPPLgrid
FRTLLIB    = $(shell gfortran.sh )


LD_FLAGS +=  -L$(LHAPDFDIR) -lLHAPDF


PDFCONVLIBS = $(shell hoppet.sh --lib  $(INSTALLBASE)/lib/libhoppet_v1.a  -L$(INSTALLBASE)/lib -lhoppet_v1 ) $(F90LIB) $(LHAPDFLIBS)

# macros to setup contact interactions and reference histograms
LHA = -DUSELHAPDF


LIBS += $(ROOTLIBS)  $(LHAPDFLIBS) $(APPLLIBS) $(FRTLLIB) $(PDFCONVLIBS)  

install : all
all : stand fnlo
stand: stand.o InitMyPdf.o
	$(CXX)  -o $@ $(LDFLAGS) $(ALLMACRO)  stand.o  InitMyPdf.o $(LIBS) 

fnlo: fnmain.o InitMyPdf.o
	$(CXX)  -o $@ $(LDFLAGS) $(ALLMACRO)  stand.o  InitMyPdf.o $(LIBS) 


.SUFFIXES : .cxx .o .f .c

.f.o :
	$(F77) $(FFLAGS) -c $<

.cxx.o:	 
	$(CXX) -c $< $(CXXFLAGS) $(ALLMACRO) -I$(INSTALLBASE)/include -I$(INSTALLBASE)/include/lhpdf $(APPLCXXFLAGS)



clean:
	rm -rf ./.libs ./.obj *.lo *.o *.la stand
