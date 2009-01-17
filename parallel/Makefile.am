AM_CPPFLAGS = -I$(top_srcdir)/include $(NETGEN_INCLUDES)

lib_LTLIBRARIES = libparallel.la

libparallel_la_SOURCES = bvpparallel.cpp parallelblockjacobi.cpp  \
paralleldofs.cpp parallelmeshaccess.cpp parallelngs.cpp	      \
parallelsparsematrix.cpp parallelvvector.cpp superlu_dist.cpp


# libparallel_la_LIBADD = $(top_builddir)/basiclinalg/libngbla.la

include_HEADERS = paralleldofs.hpp parallelmeshaccess.hpp \
parallelngs.hpp superlu_dist.hpp


libparallel_la_LDFLAGS = -avoid-version

