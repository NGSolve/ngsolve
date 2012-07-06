AM_CPPFLAGS = -I$(top_srcdir)/include $(NETGEN_INCLUDES) $(LAPACK_INCLUDES) $(PARDISO_INCLUDES) $(TCL_INCLUDES) $(MPI_INCLUDES)

lib_LTLIBRARIES = libngsolve.la

libngsolve_la_SOURCES = ngsolve.cpp bvp.cpp generalizedEVP.cpp \
numproc.cpp shapetester.cpp evplapack.cpp hyperbolic.cpp  \
nptests.cpp numprocee.cpp pde.cpp pdeparser.cpp


libngsolve_la_LIBADD = $(top_builddir)/comp/libngcomp.la	    \
$(top_builddir)/ngstd/libngstd.la 		\
$(top_builddir)/fem/libngfem.la 		\
$(top_builddir)/parallel/libparallel.la		\
$(top_builddir)/linalg/libngla.la		\
$(top_builddir)/basiclinalg/libngbla.la		\
$(top_builddir)/multigrid/libngmg.la		\
$(LAPACK_LIBS) $(TCL_LIB_SPEC)


include_HEADERS = solve.hpp pde.hpp numproc.hpp

libngsolve_la_LDFLAGS = -avoid-version -rdynamic


