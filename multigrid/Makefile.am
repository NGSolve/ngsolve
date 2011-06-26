AM_CPPFLAGS = -I$(top_srcdir)/include $(NETGEN_INCLUDES) $(MPI_INCLUDES)

lib_LTLIBRARIES = libngmg.la

libngmg_la_SOURCES = mgpre.cpp  prolongation.cpp  smoother.cpp

# libngmg_la_LIBADD = $(top_builddir)/ngstd/libngcomp.la

include_HEADERS = mgpre.hpp  multigrid.hpp  prolongation.hpp  smoother.hpp

libngmg_la_LDFLAGS = -avoid-version
