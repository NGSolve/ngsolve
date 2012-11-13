AM_CPPFLAGS = -I$(top_srcdir)/include $(MUMPS_INCLUDES) $(HYPRE_INCLUDES)
#  -I/opt/hypre-2.8.0b/include

lib_LTLIBRARIES = libparallel.la

libparallel_la_SOURCES = parallelvvector.cpp parallel_matrices.cpp

# libparallel_la_LIBADD = $(top_builddir)/comp/libngcomp.la
# libparallel_la_LIBADD = /opt/hypre-2.8.0b/lib hy/o$(top_builddir)/comp/libngcomp.la

include_HEADERS = parallelngs.hpp parallelvector.hpp parallel_matrices.hpp


libparallel_la_LDFLAGS = -avoid-version


