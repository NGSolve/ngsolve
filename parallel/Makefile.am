AM_CPPFLAGS = -I$(top_srcdir)/include $(NETGEN_INCLUDES) $(MPI_INCLUDES)  $(MUMPS_INCLUDES)

lib_LTLIBRARIES = libparallel.la

libparallel_la_SOURCES =   paralleldofs.cpp  \
parallelvvector.cpp parallel_matrices.cpp

# libparallel_la_LIBADD = $(top_builddir)/comp/libngcomp.la

include_HEADERS = paralleldofs.hpp parallelngs.hpp \
parallelvector.hpp parallel_matrices.hpp mpiwrapper.hpp


libparallel_la_LDFLAGS = -avoid-version

