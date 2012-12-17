AM_CPPFLAGS = -I$(top_srcdir)/include $(PARDISO_INCLUDES) $(TCL_INCLUDES)

lib_LTLIBRARIES = libngsolve.la
bin_PROGRAMS = ngs




libngsolve_la_SOURCES = bvp.cpp generalizedEVP.cpp \
numproc.cpp evplapack.cpp hyperbolic.cpp  \
nptests.cpp numprocee.cpp pde.cpp pdeparser.cpp

if NGSGUI
  libngsolve_la_SOURCES += ngsolve.cpp shapetester.cpp
  AM_CPPFLAGS += -DTCL
endif

libngsolve_la_LIBADD = $(top_builddir)/comp/libngcomp.la	    \
$(top_builddir)/ngstd/libngstd.la 		\
$(top_builddir)/fem/libngfem.la 		\
$(top_builddir)/parallel/libparallel.la		\
$(top_builddir)/linalg/libngla.la		\
$(top_builddir)/basiclinalg/libngbla.la		\
$(top_builddir)/multigrid/libngmg.la		\
$(LAPACK_LIBS) $(TCL_LIB_SPEC)

libngsolve_la_LDFLAGS = -avoid-version -rdynamic



ngs_SOURCES = ngs.cpp
ngs_LDADD = libngsolve.la \
	$(top_builddir)/comp/libngcomp.la \
	$(top_builddir)/multigrid/libngmg.la \
	$(top_builddir)/fem/libngfem.la \
	$(top_builddir)/parallel/libparallel.la \
	$(top_builddir)/linalg/libngla.la \
	$(top_builddir)/ngstd/libngstd.la \
	$(LAPACK_LIBS) \
        -L$(libdir) -lnglib -linterface 



include_HEADERS = solve.hpp pde.hpp numproc.hpp



