AM_CPPFLAGS = -I$(top_builddir) -I$(top_srcdir)/include -I/usr/include/python$(PYTHON_VERSION)

lib_LTLIBRARIES = libngbla.la

libngbla_la_SOURCES = bandmatrix.cpp calcinverse.cpp cholesky.cpp \
eigensystem.cpp vecmat.cpp LapackGEP.cpp LapackInterface.hpp python_bla.cpp 


libngbla_la_LIBADD = $(top_builddir)/ngstd/libngstd.la $(LAPACK_LIBS) $(PYTHON_LIBS)
# libngbla_la_LIBADD += -lgmp -lgmpxx

include_HEADERS = bandmatrix.hpp cholesky.hpp matrix.hpp ng_lapack.hpp \
vector.hpp bla.hpp expr.hpp symmetricmatrix.hpp arch.hpp clapack.h \
tensor.hpp cuda_bla.hpp

libngbla_la_LDFLAGS = -avoid-version 

