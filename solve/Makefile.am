AM_CPPFLAGS = -I$(top_srcdir)/include $(NETGEN_INCLUDES) $(LAPACK_INCLUDES) $(PARDISO_INCLUDES) $(TCL_INCLUDES)
# -I/home/joachim/netgen-mesher/netgen/libsrc/interface

lib_LTLIBRARIES = libngsolve.la

libngsolve_la_SOURCES = ngsolve.cpp bvp.cpp generalizedEVP.cpp \
numproc.cpp pde.cpp shapetester.cpp evplapack.cpp hyperbolic.cpp  \
nptests.cpp numprocee.cpp pdeparser.cpp


libngsolve_la_LIBADD = $(top_builddir)/comp/libngcomp.la	    \
$(top_builddir)/fem/libngfemng.la $(top_builddir)/linalg/libngla.la \
$(top_builddir)/basiclinalg/libngbla.la				    \
$(top_builddir)/multigrid/libngmg.la  $(PARDISO_LIBS) $(MUMPS_LIBS) $(LAPACK_LIBS) 


include_HEADERS = numproc.hpp pde.hpp solve.hpp 

libngsolve_la_LDFLAGS = -rdynamic
# libngsolve_la_LDFLAGS = -avoid-version -rdynamic


