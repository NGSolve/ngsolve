#
#  libngfemng ... the library linked with Netgen/NGSolve
#  libngfem ..... a stand-alone fem library
#


lib_LTLIBRARIES = libngfemng.la libngfem.la

libngfemng_la_CPPFLAGS = -I$(top_srcdir)/include  $(NETGEN_INCLUDES) $(LAPACK_INCLUDES)

libngfemng_la_LDFLAGS = -avoid-version -rdynamic 

libngfemng_la_SOURCES = bdbequations.cpp recursive_pol_trig.cpp	       \
coefficient.cpp integrator.cpp specialelement.cpp elementtopology.cpp  \
hcurlfe.cpp intrule.cpp facetfe.cpp hcurlhofe.cpp l2hofe.cpp	       \
fastmat.cpp vectorfacetfe.cpp finiteelement.cpp scalarfe.cpp	       \
hdiv_equations.cpp maxwellintegrator.cpp generic_recpol.cpp hdivfe.cpp \
h1hofe.cpp h1hofefo.cpp hdivhofe.cpp recursive_pol.cpp hybridDG.cpp    \
pml.cpp


libngfem_la_CPPFLAGS = -I$(top_srcdir)/include  $(LAPACK_INCLUDE)
libngfem_la_LDFLAGS = -avoid-version

libngfem_la_SOURCES = bdbequations.cpp recursive_pol_trig.cpp	       \
coefficient.cpp integrator.cpp specialelement.cpp elementtopology.cpp  \
hcurlfe.cpp intrule.cpp facetfe.cpp hcurlhofe.cpp l2hofe.cpp	       \
fastmat.cpp vectorfacetfe.cpp finiteelement.cpp scalarfe.cpp	       \
hdiv_equations.cpp maxwellintegrator.cpp generic_recpol.cpp hdivfe.cpp \
h1hofe.cpp h1hofefo.cpp hdivhofe.cpp recursive_pol.cpp dynamicmem.cpp  \
hybridDG.cpp


include_HEADERS = bdbequations.hpp bdbintegrator.hpp coefficient.hpp   \
elementtopology.hpp elementtransformation.hpp facetfe.hpp fem.hpp      \
femtypes.hpp finiteelement.hpp scalarfe.hpp generic_recpol.hpp	       \
h1hofe.hpp h1hofefo.hpp hcurlfe.hpp hcurlhofe.hpp hdiv_equations.hpp   \
hdivfe.hpp hdivhofe.hpp integrator.hpp intrule.hpp l2hofe.hpp	       \
recursive_pol.hpp recursive_pol_tet.hpp recursive_pol_trig.hpp	       \
 specialelement.hpp vectorfacetfe.hpp pml.hpp			       