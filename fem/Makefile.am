AM_CPPFLAGS = -I$(top_srcdir)/include  $(NETGEN_INCLUDES) $(LAPACK_INCLUDE)

noinst_LTLIBRARIES = libngfemng.la
libngfemng_la_LIBADD = $(top_builddir)/basiclinalg/libngbla.la  $(LAPACK_LIB)
libngfemng_la_LDFLAGS = -avoid-version

libngfemng_la_SOURCES = bdbequations.cpp recursive_pol_trig.cpp	      \
coefficient.cpp h1hofetp.cpp integrator.cpp specialelement.cpp	      \
elementtopology.cpp hcurlfe.cpp intrule.cpp facetfe.cpp hcurlhofe.cpp \
l2hofe.cpp fastmat.cpp hcurlhofetp.cpp l2hofetp.cpp vectorfacetfe.cpp \
finiteelement.cpp hdiv_equations.cpp maxwellintegrator.cpp	      \
generic_recpol.cpp hdivfe.cpp pml.cpp h1hofe.cpp hdivhofe.cpp	      \
recursive_pol.cpp





lib_LTLIBRARIES = libngfem.la
libngfem_la_CPPFLAGS = -I$(top_srcdir)/include  $(LAPACK_INCLUDE)
libngfem_la_LDFLAGS = -avoid-version  $(LAPACK_LIB)

libngfem_la_SOURCES = bdbequations.cpp recursive_pol_trig.cpp	      \
coefficient.cpp h1hofetp.cpp integrator.cpp specialelement.cpp	      \
elementtopology.cpp hcurlfe.cpp intrule.cpp facetfe.cpp hcurlhofe.cpp \
l2hofe.cpp fastmat.cpp hcurlhofetp.cpp l2hofetp.cpp vectorfacetfe.cpp \
finiteelement.cpp hdiv_equations.cpp maxwellintegrator.cpp	      \
generic_recpol.cpp hdivfe.cpp h1hofe.cpp hdivhofe.cpp		      \
recursive_pol.cpp dynamicmem.cpp



include_HEADERS = bdbequations.hpp bdbintegrator.hpp coefficient.hpp \
elementtopology.hpp elementtransformation.hpp facetfe.hpp fem.hpp    \
femtypes.hpp finiteelement.hpp generic_recpol.hpp h1hofe.hpp	     \
h1hofetp.hpp hcurlfe.hpp hcurlhofe.hpp hcurlhofetp.hpp		     \
hdiv_equations.hpp hdivfe.hpp hdivhofe.hpp highorderfe.hpp	     \
integrator.hpp intrule.hpp l2hofe.hpp l2hofetp.hpp pml.hpp	     \
recursive_pol.hpp recursive_pol_tet.hpp recursive_pol_trig.hpp	     \
specialelement.hpp vectorfacetfe.hpp

