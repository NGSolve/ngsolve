AM_CPPFLAGS = -I$(top_srcdir)/include  $(NETGEN_INCLUDES)

lib_LTLIBRARIES = libngfem.la

libngfem_la_SOURCES = bdbequations.cpp recursive_pol_trig.cpp	       \
coefficient.cpp h1hofetp.cpp integrator.cpp specialelement.cpp	       \
elementtopology.cpp hcurlfe.cpp intrule.cpp tpequations.cpp	       \
facetfe.cpp hcurlhofe.cpp l2hofe.cpp tpmaxwell.cpp fastmat.cpp	       \
hcurlhofetp.cpp l2hofetp.cpp vectorfacetfe.cpp finiteelement.cpp       \
hdiv_equations.cpp maxwellintegrator.cpp generic_recpol.cpp hdivfe.cpp \
pml.cpp h1hofe.cpp hdivhofe.cpp recursive_pol.cpp


libngfem_la_LIBADD = $(top_builddir)/basiclinalg/libngbla.la

include_HEADERS = bdbequations.hpp bdbintegrator.hpp coefficient.hpp elementtopology.hpp elementtransformation.hpp facetfe.hpp fem.hpp femtypes.hpp finiteelement.hpp generic_recpol.hpp h1hofe.hpp h1hofetp.hpp h1hofetp_save.hpp hcurlfe.hpp hcurlhofe.hpp hcurlhofetp.hpp hdiv_equations.hpp hdivfe.hpp hdivhofe.hpp highorderfe.hpp integrator.hpp intrule.hpp l2hofe.hpp l2hofetp.hpp pml.hpp recursive_pol.hpp recursive_pol_tet.hpp recursive_pol_trig.hpp specialelement.hpp vectorfacetfe.hpp

libngfem_la_LDFLAGS = -avoid-version

