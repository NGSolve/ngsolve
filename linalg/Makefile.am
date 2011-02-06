AM_CPPFLAGS = -I$(top_srcdir)/include $(NETGEN_INCLUDES) $(LAPACK_INCLUDES) $(PARDISO_INCLUDES) $(MUMPS_INCLUDES) $(MPI_INCLUDES)

lib_LTLIBRARIES = libngla.la

libngla_la_SOURCES = basematrix.cpp basevector.cpp blockjacobi.cpp cg.cpp chebyshev.cpp commutingAMG.cpp eigen.cpp jacobi.cpp order.cpp pardisoinverse.cpp sparsecholesky.cpp sparsematrix.cpp special_matrix.cpp superluinverse.cpp mumpsinverse.cpp vvector.cpp elementbyelement.cpp

libngla_la_LIBADD = $(top_builddir)/basiclinalg/libngbla.la 
# $(LAPACK_LIB)

include_HEADERS = basematrix.hpp basevector.hpp blockjacobi.hpp cg.hpp chebyshev.hpp commutingAMG.hpp eigen.hpp jacobi.hpp la.hpp order.hpp pardisoinverse.hpp sparsecholesky.hpp sparsematrix.hpp special_matrix.hpp superluinverse.hpp mumpsinverse.hpp vvector.hpp elementbyelement.hpp

libngla_la_LDFLAGS = -avoid-version $(MUMPS_LIBS)
