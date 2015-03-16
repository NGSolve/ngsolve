AM_CPPFLAGS = -I$(top_builddir) $(GOLD_INCLUDE) -I$(srcdir) -I$(top_srcdir)/include 

lib_LTLIBRARIES = libngfem.la

libngfem_la_LIBADD = $(top_builddir)/ngstd/libngstd.la \
	$(top_builddir)/basiclinalg/libngbla.la $(LAPACK_LIBS)

libngfem_la_SOURCES = bdbequations.cpp maxwellintegrator.cpp  \
hdiv_equations.cpp h1hofe.cpp h1lofe.cpp l2hofe.cpp hcurlhofe.cpp     \
hcurlfe.cpp vectorfacetfe.cpp hdivhofe.cpp recursive_pol_trig.cpp     \
coefficient.cpp integrator.cpp specialelement.cpp elementtopology.cpp \
intrule.cpp fastmat.cpp finiteelement.cpp elementtransformation.cpp   \
scalarfe.cpp generic_recpol.cpp hdivfe.cpp recursive_pol.cpp	      \
hybridDG.cpp diffop.cpp l2hofefo.cpp h1hofefo.cpp   \
facethofe.cpp python_fem.cpp fem_kernels.cu  DGIntegrators.cpp pml.cpp \
h1hofe_segm.cpp h1hofe_trig.cpp
# 
#  

include_HEADERS = bdbequations.hpp hcurl_equations.hpp		       \
elasticity_equations.hpp diffop.hpp bdbintegrator.hpp coefficient.hpp  \
elementtopology.hpp elementtransformation.hpp facetfe.hpp fastmat.hpp  \
facethofe.hpp fem.hpp finiteelement.hpp generic_recpol.hpp	       \
h1hofefo.hpp h1hofefo_impl.hpp h1hofe.hpp h1lofe.hpp hcurlfe.hpp       \
hcurlhofe.hpp hcurllofe.hpp hdiv_equations.hpp hdivfe.hpp hdivhofe.hpp \
integrator.hpp intrule.hpp l2hofefo.hpp l2hofe.hpp recursive_pol.hpp   \
recursive_pol_tet.hpp recursive_pol_trig.hpp scalarfe.hpp	       \
specialelement.hpp thdivfe.hpp tscalarfe.hpp vectorfacetfe.hpp	       \
hdivlofe.hpp hdivhofefo.hpp pml.hpp precomp.hpp h1hofe_impl.hpp	       \
hdivhofe_impl.hpp tscalarfe_impl.hpp thdivfe_impl.hpp l2hofe_impl.hpp  \
diffop_impl.hpp hcurlhofe_impl.hpp thcurlfe.hpp thcurlfe_impl.hpp


libngfem_la_LDFLAGS = -avoid-version



if NGSCUDA


# dost not yet work. we need the legendre coefficients extern, so we need -dc cuda compilation
# 
NVCC = /usr/local/cuda-6.5/bin/nvcc -ccbin=gcc-4.8
SUFFIXES = .cu
.cu.o:
	$(NVCC) -arch compute_35 -code sm_35 -std=c++11 -O2 -I$(top_srcdir)/include  -rdc true -Xcompiler -fpic -DUSE_MYCOMPLEX -DHD="__host__ __device__" -c -o tmp.o $<
	$(NVCC) -arch compute_35 -code sm_35 -Xcompiler -fpic -dlink tmp.o -o $@ -L/usr/local/cuda-6.5/lib64 -lcudadevrt -lcudart

$(top_builddir)/fem/tmp.o :  $(top_builddir)/fem/fem_kernels.o

libngfem_la_LIBADD += $(top_builddir)/fem/fem_kernels.o $(top_builddir)/fem/tmp.o -L/usr/local/cuda-6.5/lib64  -lcudadevrt -lcudart

endif

dist_bin_SCRIPTS = ngfemcf.py
