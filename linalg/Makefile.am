
AM_CPPFLAGS = -I$(top_builddir)  -I$(top_srcdir)/include $(PARDISO_INCLUDES) $(MUMPS_INCLUDES) $(SUPERLU_INCLUDES)

lib_LTLIBRARIES = libngla.la

libngla_la_SOURCES = linalg_kernels.cu basematrix.cpp basevector.cpp \
blockjacobi.cpp cg.cpp chebyshev.cpp commutingAMG.cpp eigen.cpp	     \
jacobi.cpp order.cpp pardisoinverse.cpp sparsecholesky.cpp	     \
sparsematrix.cpp special_matrix.cpp superluinverse.cpp		     \
mumpsinverse.cpp elementbyelement.cpp arnoldi.cpp paralleldofs.cpp   \
cuda_linalg.cpp python_linalg.cpp

libngla_la_LIBADD = $(top_builddir)/basiclinalg/libngbla.la \
  $(top_builddir)/ngstd/libngstd.la \
  $(top_builddir)/parallel/libparallel.la \
  $(LAPACK_LIBS) 
#  -L$(libdir) -lnglib -linterface


if NGSCUDA

NVCC = /usr/local/cuda-6.5/bin/nvcc -ccbin=gcc-4.8
SUFFIXES = .cu
.cu.o:
	$(NVCC) -Xcompiler -fpic -c -o $@ $<

libngla_la_LIBADD += $(top_builddir)/linalg/linalg_kernels.o

endif



include_HEADERS = basematrix.hpp basevector.hpp blockjacobi.hpp cg.hpp \
chebyshev.hpp commutingAMG.hpp eigen.hpp jacobi.hpp la.hpp order.hpp   \
pardisoinverse.hpp sparsecholesky.hpp sparsematrix.hpp		       \
special_matrix.hpp superluinverse.hpp mumpsinverse.hpp vvector.hpp     \
elementbyelement.hpp arnoldi.hpp paralleldofs.hpp cuda_linalg.hpp  \
umfpackinverse.hpp

libngla_la_LDFLAGS = -avoid-version $(PARDISO_LIBS) $(MUMPS_LIBS) \
$(SUPERLU_LIBS) $(LAPACK_LIBS) $(PYTHON_LIBS)


# dist_bin_SCRIPTS = la.py
