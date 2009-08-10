AUTOMAKE_OPTIONS = foreign
LHA_VERSION = 1:0:0

#  Header files
includedir= @includedir@/lhpdf

include_HEADERS = lhpdf.h \
\
lhpdf_alekhin.h lhpdf_botje.h lhpdf_cteq6.h \
lhpdf_fermi2002.h lhpdf_mrst2001.h

#  Core library : liblhpdf.la
lib_LTLIBRARIES = liblhpdf.la liblhpdf_cteq6.la liblhpdf_mrst2001.la \
liblhpdf_fermi2002.la liblhpdf_alekhin.la liblhpdf_botje.la

liblhpdf_la_SOURCES = lhpdf.c lhpdfset.c evlcteq.c \
f77-evlcteq.f Sqcdnum.f f77-qcdnum.f evlcteq-c.f qcdnum-c.f qcdnum.c

liblhpdf_la_LDFLAGS = -version-info $(LHA_VERSION)
liblhpdf_la_LIBADD = $(FLIBS) 


#    Libraries for the different pdfsets
liblhpdf_cteq6_la_SOURCES = cteq6.c
liblhpdf_cteq6_la_LDFLAGS = -version-info $(LHA_VERSION)
liblhpdf_cteq6_la_LIBADD = liblhpdf.la 

liblhpdf_mrst2001_la_SOURCES = mrst2001.c
liblhpdf_mrst2001_la_LDFLAGS = -version-info $(LHA_VERSION)
liblhpdf_mrst2001_la_LIBADD = liblhpdf.la 

liblhpdf_fermi2002_la_SOURCES = fermi2002.c 
liblhpdf_fermi2002_la_LDFLAGS = -version-info $(LHA_VERSION)
liblhpdf_fermi2002_la_LIBADD = liblhpdf.la 

liblhpdf_alekhin_la_SOURCES = alekhin.c 
liblhpdf_alekhin_la_LDFLAGS = -version-info $(LHA_VERSION)
liblhpdf_alekhin_la_LIBADD = liblhpdf.la 

liblhpdf_botje_la_SOURCES = botje.c
liblhpdf_botje_la_LDFLAGS = -version-info $(LHA_VERSION)
liblhpdf_botje_la_LIBADD = liblhpdf.la 


