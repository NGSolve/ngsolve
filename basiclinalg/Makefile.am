AM_CPPFLAGS = -I$(top_srcdir)/include $(NETGEN_INCLUDES) $(LAPACK_INCLUDES)

lib_LTLIBRARIES = libngbla.la

libngbla_la_SOURCES = bandmatrix.cpp calcinverse.cpp cholesky.cpp \
eigensystem.cpp vecmat.cpp LapackGEP.cpp LapackInterface.hpp

libngbla_la_LIBADD = $(top_builddir)/ngstd/libngstd.la

include_HEADERS = bandmatrix.hpp cholesky.hpp matrix.hpp ng_lapack.hpp \
vector.hpp bla.hpp expr.hpp symmetricmatrix.hpp

libngbla_la_LDFLAGS = -avoid-version 
