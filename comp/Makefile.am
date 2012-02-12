AM_CPPFLAGS = -I$(top_srcdir)/include $(NETGEN_INCLUDES) $(LAPACK_INCLUDES) $(TCL_INCLUDES) $(MPI_INCLUDES)

lib_LTLIBRARIES = libngcomp.la

libngcomp_la_SOURCES = bilinearform.cpp facetfespace.cpp fespace.cpp \
gridfunction.cpp h1hofespace.cpp hcurlhdivfes.cpp hcurlhofespace.cpp \
hdivfes.cpp hdivhofespace.cpp hierarchicalee.cpp l2hofespace.cpp     \
linearform.cpp meshaccess.cpp ngsobject.cpp postproc.cpp	     \
preconditioner.cpp vectorfacetfespace.cpp bddc.cpp 

libngcomp_la_LIBADD = $(top_builddir)/fem/libngfem.la \
$(top_builddir)/linalg/libngla.la \
$(top_builddir)/parallel/libparallel.la \
$(top_builddir)/multigrid/libngmg.la	\
$(top_builddir)/ngstd/libngstd.la		

include_HEADERS = bilinearform.hpp comp.hpp facetfespace.hpp	       \
fespace.hpp gridfunction.hpp h1hofespace.hpp hcurlhdivfes.hpp	       \
hcurlhofespace.hpp hdivfes.hpp hdivhofespace.hpp		       \
 l2hofespace.hpp linearform.hpp meshaccess.hpp ngsobject.hpp	       \
 postproc.hpp preconditioner.hpp vectorfacetfespace.hpp wirebasket.hpp \
bddc.hpp 

libngcomp_la_LDFLAGS = -avoid-version
