AM_CPPFLAGS = -I$(top_srcdir)/include $(NETGEN_INCLUDES)

lib_LTLIBRARIES = libngcomp.la

libngcomp_la_SOURCES = bilinearform.cpp facetfespace.cpp fespace.cpp \
gridfunction.cpp h1hofespace.cpp hcurlhdivfes.cpp hcurlhofespace.cpp \
hdivfes.cpp hdivhofespace.cpp hdivhybridhofespace.cpp		     \
hierarchicalee.cpp highorderfes.cpp l2hofespace.cpp linearform.cpp   \
meshaccess.cpp ngsobject.cpp postproc.cpp preconditioner.cpp	     \
 vectorfacetfespace.cpp

# libngcomp_la_LIBADD = $(top_builddir)/ngstd/libngfem.la

include_HEADERS = bilinearform.hpp blocked_gridfunction.hpp comp.hpp \
facetfespace.hpp fespace.hpp gridfunction.hpp h1hofespace.hpp	     \
hcurlhdivfes.hpp hcurlhofespace.hpp hdivfes.hpp hdivhofespace.hpp    \
hdivhybridhofespace.hpp highorderfes.hpp l2hofespace.hpp	     \
linearform.hpp meshaccess.hpp ngsobject.hpp postproc.hpp	     \
preconditioner.hpp vectorfacetfespace.hpp


libngcomp_la_LDFLAGS = -avoid-version
