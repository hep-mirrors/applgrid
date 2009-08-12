AUTOMAKE_OPTIONS = foreign
SUBDIRS = libltdl nlo-hep nlo-core nlo-ps nlo-amp proc-epa proc-dis proc-hhc \
proc-hhc2ph proc-photo nlo-phys . src

SUBLIBS = nlo-hep/hep-lorentz/libhep-lorentz.la nlo-hep/hep-random/libhep-random.la \
nlo-hep/hep-stat/libhep-stat.la \
\
nlo-core/libcore-nlo++.la \
nlo-ps/libps-nlo++.la nlo-phys/libphys-nlo++.la nlo-amp/libamp-nlo++.la \
proc-epa/libepa-nlo++.la proc-dis/libdis-nlo++.la \
proc-hhc/libhhc-nlo++.la proc-hhc2ph/libhhc2ph-nlo++.la \
proc-photo/libphoto-nlo++.la



lib_LTLIBRARIES = libnlo++.la
libnlo___la_SOURCES = version.cc
libnlo___la_LDFLAGS = -version-info $(LIB_VERSION)
libnlo___la_LIBADD = $(SUBLIBS)



libtool-nlojet++: libtool
	cp libtool libtool-nlojet++


install-exec-local: libtool-nlojet++
	$(mkinstalldirs) $(DESTDIR)$(libexecdir)/nlo++
	$(INSTALL_SCRIPT) libtool-nlojet++ $(DESTDIR)$(libexecdir)/nlo++



