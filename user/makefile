
CXX = g++
F77 = gfortran

FFLAGS += -g 


CXXFLAGS +=  -O3 -fPIC -g 

# root
ROOTINCS = $(shell root-config --cflags) 
ROOTLIBS = $(shell root-config --glibs) 
ROOTARCH = $(findstring -m64, $(ROOTINCS) )


#LHAPDF
LHAPDFINCS = -I$(shell lhapdf-config --prefix)/include
LHAPDFDIR  = $(shell lhapdf-config --prefix)/lib
LHAPDFLIBS = -L$(LHAPDFDIR) -lLHAPDF

# applgrid
APPLCXXFLAGS = $(shell applgrid-config --cxxflags)
APPLCLIBS    = $(shell applgrid-config --ldcflags)
APPLLIBS    = $(shell applgrid-config --ldflags)

# hoppet
HOPPETLIBS =  $(shell hoppet-config --libs)

FRTLIB   = $(shell gfortran $(CXXFLAGS) -print-file-name=libgfortran.a)
FRTLLIB  = -L$(subst /libgfortran.a, ,$(FRTLIB) ) -lgfortran 

LD_FLAGS += $(ROOTARCH) -L$(LHAPDFDIR) -lLHAPDF

CXXFLAGS += $(ROOTARCH) $(ROOTINCS) $(APPLCXXFLAGS) $(LHAPDFINCS)
LDFLAGS  += $(ROOTARCH)
FFLAGS   += $(ROOTARCH)


FRTLLIB = $(shell gfortran $(CXXFLAGS) -print-file-name=libgfortran.a)
FRTLIB  = -L$(subst /libgfortran.a, ,$(FRTLLIB) ) -lgfortran

LD_FLAGS += $(ROOTARCH) $(LHAPDFLIBS) $(APPLLIBS) $(HOPPETLIBS) $(FRTLIB)


PDFCONVLIBS = $(HOPPETLIBS) $(LHAPDFLIBS)

# macros to setup contact interactions and reference histograms
LHA = -DUSELHAPDF


LIBS += $(ROOTLIBS)  $(LHAPDFLIBS) $(FRTLIB)  $(FRTLIB) $(PDFCONVLIBS)  

install : all
all : stand fnlo fstand

stand: stand.o 
	$(CXX)  -o $@ $(LDFLAGS) $(ALLMACRO) $<  $(LIBS) $(APPLCLIBS) 


fstand: fstand.o fmain.o
	@ echo $(LIBS)
	$(CXX)  -o $@ fstand.o fmain.o  $(LDFLAGS) $(ALLMACRO)  $(LIBS) $(APPLLIBS) 


fnlo: fnmain.o 
	$(CXX)  -o $@ $(LDFLAGS) $(ALLMACRO) $<  $(LIBS) $(APPLCLIBS) 


.SUFFIXES : .cxx .o .f .c

.f.o :
	$(F77) $(FFLAGS) $(CXXFLAGS) -c $<

.cxx.o:	 
	$(CXX) -c $< $(CXXFLAGS) $(ALLMACRO)  $(APPLCXXFLAGS)



clean:
	rm -rf ./.libs ./.obj *.lo *.o *.la stand fnlo fstand *~
