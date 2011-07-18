AM_CPPFLAGS = -I$(top_srcdir)/include $(NETGEN_INCLUDES) $(MPI_INCLUDES)  $(MUMPS_INCLUDES)

lib_LTLIBRARIES = libparallel.la

libparallel_la_SOURCES =   paralleldofs.cpp  \
parallelngs.cpp parallelvvector.cpp parallel_matrices.cpp

# parallelblockjacobi.cpp parallelsparsematrix.cpp superlu_dist.cpp parallelmeshaccess.cpp

libparallel_la_LIBADD = $(top_builddir)/comp/libngcomp.la

include_HEADERS = paralleldofs.hpp parallelngs.hpp \
parallelvector.hpp parallel_matrices.hpp mpiwrapper.hpp

#  superlu_dist.hpp parallelsparsematrix.hpp parallelmeshaccess.hpp

libparallel_la_LDFLAGS = -avoid-version

