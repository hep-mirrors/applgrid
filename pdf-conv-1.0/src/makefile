# Makefile generated automatically with
# /home/salam/scripts/makePNEW.perl libdglap.a
# default program to compile
PROG =	libdglap.a

ALLPROG = 	

# INSTALLDIR = $(HOME)

ifdef BASEDIR
INSTALLBASE = $(BASEDIR)
else
INSTALLBASE = ../..
endif


# This will be used one day...
ALLPROGSRC =	

ALLPROGOBJ =	

SRCS =	hplog.f xpij2e.f xpns2e.f as_server.f90 assertions.f90 \
	coefficient_functions.f90 conv_objects.f90 convolution.f90 \
	dglap_choices.f90 evolution.f90 f77_pdftab.f90 holders.f90 \
	integrator.f90 interpolation.f90 new_as.f90 pdf_general.f90 \
	pdf_representation.f90 pdf_tabulate.f90 qcd.f90 runge_kutta.f90 \
	sort.f90 special_functions.f90 splitting_functions.f90 \
	splitting_functions_nnlo.f90 types.f90 warnings_and_errors.f90 \
	xa2hgp.f90 xpij2n.f90 xpij2p.f90 xpns2n.f90 xpns2p.f90 \
	dglapInit.f dglapEvolve.f

POSTSRCS =	

OBJS =	hplog.o xpij2e.o xpns2e.o as_server.o assertions.o \
	coefficient_functions.o conv_objects.o convolution.o dglap_choices.o \
	evolution.o f77_pdftab.o holders.o integrator.o interpolation.o \
	new_as.o pdf_general.o pdf_representation.o pdf_tabulate.o qcd.o \
	runge_kutta.o sort.o special_functions.o splitting_functions.o \
	splitting_functions_nnlo.o types.o warnings_and_errors.o xa2hgp.o \
	xpij2n.o xpij2p.o xpns2n.o xpns2p.o dglapInit.o dglapEvolve.o

POSTOBJS =	
POSTLIB = 

LIBS =  

CC = cc
CFLAGS += -O3 
#FC = lf95
F77 = gfortran
#F90 = lf95
F90 = gfortran
# F90FLAGS = -m32


AR = ar
ARFLAGS = -r

RANLIB = ranlib


CP = cp


all: libdglap.a

libdglap.a: $(OBJS)
	$(AR) $(ARFLAGS) $@ $(OBJS)
	$(RANLIB) $@

ALL: libdglap.a  $(ALLPROG)

libclean:
	rm -f  $(ALLPROGOBJS) $(OBJS) $(POSTOBJS) 

clean:
	rm -f  $(ALLPROGOBJS) $(OBJS) $(POSTOBJS) *.mod *.d *.a

realclean:
	rm -f libdglap.a $(ALLPROG) $(ALLPROGOBJ) $(OBJS) $(POSTOBJS) *.mod *.d


install: libdglap.a
	./mkdir.sh $(INSTALLDIR)/lib
	$(CP) $< $(INSTALLDIR)/lib


.SUFFIXES: $(SUFFIXES) .f90

%.o: %.f90 
	$(F90) $(F90FLAGS) -c $<

%.o: %.f 
	$(F77) $(FFLAGS) -c $<

xpij2e.o: qcd.o

xpns2e.o: qcd.o 

as_server.o: assertions.o new_as.o types.o warnings_and_errors.o 

assertions.o: types.o

coefficient_functions.o: convolution.o qcd.o types.o

conv_objects.o: assertions.o convolution.o dglap_choices.o \
	pdf_representation.o qcd.o splitting_functions.o types.o \
	warnings_and_errors.o

convolution.o: assertions.o integrator.o interpolation.o sort.o types.o \
	warnings_and_errors.o

evolution.o: as_server.o assertions.o conv_objects.o convolution.o \
	dglap_choices.o holders.o pdf_representation.o qcd.o runge_kutta.o \
	types.o warnings_and_errors.o

f77_pdftab.o: conv_objects.o convolution.o dglap_choices.o holders.o \
	pdf_general.o pdf_tabulate.o types.o warnings_and_errors.o

holders.o: assertions.o coefficient_functions.o conv_objects.o convolution.o \
	dglap_choices.o pdf_representation.o qcd.o types.o \
	warnings_and_errors.o

integrator.o: types.o

interpolation.o: assertions.o types.o warnings_and_errors.o

new_as.o: assertions.o qcd.o runge_kutta.o types.o warnings_and_errors.o

pdf_general.o: convolution.o pdf_representation.o types.o

pdf_representation.o: assertions.o types.o warnings_and_errors.o

pdf_tabulate.o: as_server.o conv_objects.o convolution.o evolution.o \
	holders.o interpolation.o pdf_general.o pdf_representation.o types.o \
	warnings_and_errors.o

qcd.o: types.o

runge_kutta.o: types.o

sort.o: assertions.o types.o warnings_and_errors.o

special_functions.o: types.o

splitting_functions.o: coefficient_functions.o convolution.o qcd.o \
	special_functions.o splitting_functions_nnlo.o types.o \
	warnings_and_errors.o

splitting_functions_nnlo.o: convolution.o dglap_choices.o qcd.o types.o \
	warnings_and_errors.o xpij2e.o xpij2n.o xpij2p.o xpns2e.o xpns2n.o \
	xpns2p.o

types.o: 

warnings_and_errors.o: types.o

archive:
	cd ../.. ; tar -czf pdf-conv-1.0.tgz pdf-conv-1.0/src/{makefile,*.f*,*.mod,*.sh}

