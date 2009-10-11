AUTOMAKE_OPTIONS = foreign
includedir= @includedir@/nlo++/bits
bindir=$(libexecdir)/nlo++


include_HEADERS = nlo-event.h nlo-innerprod.h nlo-phasespace.h \
\
nlo-split.h nlo-weight.h nlo-dipole_mom.h \
nlo-color.h nlo-flkern.h \
\
nlo-process.h nlo-process_i0f0.h nlo-process_i1f0.h nlo-process_i2f0.h \
\
nlo-jetfunc.h nlo-basic_user.h nlo-basic_user_set.h  nlo-basic_user_result.h \
\
nlo-integral.h nlo-integral_i0f0.h nlo-integral_i1f0.h nlo-integral_i2f0.h 

bin_PROGRAMS = mk-basic_user_set-header mk-basic_user_result-header
mk_basic_user_set_header_SOURCES = mkbuset.cc
mk_basic_user_result_header_SOURCES = mkburesult.cc

nlo-basic_user_set.h: mk-basic_user_set-header
	test -f $@ ||  rm -f nlo-basic_user_set.h
	./mk-basic_user_set-header > nlo-basic_user_set.h

nlo-basic_user_result.h: mk-basic_user_result-header
	test -f $@ ||  rm -f nlo-basic_user_result.h
	./mk-basic_user_result-header > nlo-basic_user_result.h


