AM_CPPFLAGS = -I$(top_srcdir)/include

lib_LTLIBRARIES = libngfem.la

libngfem_la_LDFLAGS = -avoid-version
libngfem_la_LIBADD = $(top_builddir)/ngstd/libngstd.la \
	$(top_builddir)/basiclinalg/libngbla.la $(LAPACK_LIBS)

libngfem_la_SOURCES = bdbequations.cpp pml.cpp maxwellintegrator.cpp \
hdiv_equations.cpp h1hofe.cpp h1lofe.cpp l2hofe.cpp hcurlhofe.cpp    \
hcurlfe.cpp vectorfacetfe.cpp hdivhofe.cpp \
recursive_pol_trig.cpp coefficient.cpp integrator.cpp		     \
specialelement.cpp elementtopology.cpp intrule.cpp fastmat.cpp	     \
finiteelement.cpp elementtransformation.cpp scalarfe.cpp	     \
generic_recpol.cpp hdivfe.cpp recursive_pol.cpp hybridDG.cpp	     \
DGIntegrators.cpp diffop.cpp

# l2hofefo.cpp h1hofefo.cpp

include_HEADERS = bdbequations.hpp hcurl_equations.hpp		      \
elasticity_equations.hpp diffop.hpp bdbintegrator.hpp coefficient.hpp \
elementtopology.hpp elementtransformation.hpp facetfe.hpp	      \
facethofe.hpp fem.hpp finiteelement.hpp generic_recpol.hpp	      \
h1hofefo.hpp h1hofe.hpp h1lofe.hpp hcurlfe.hpp hcurlhofe.hpp	      \
hcurllofe.hpp hdiv_equations.hpp hdivfe.hpp hdivhofe.hpp	      \
integrator.hpp intrule.hpp l2hofefo.hpp l2hofe.hpp recursive_pol.hpp  \
recursive_pol_tet.hpp recursive_pol_trig.hpp scalarfe.hpp	      \
specialelement.hpp thdivfe.hpp tscalarfe.hpp vectorfacetfe.hpp	      \
hdivlofe.hpp hdivhofefo.hpp pml.hpp precomp.hpp h1hofe_impl.hpp	      \
hdivhofe_impl.hpp tscalarfe_impl.hpp thdivfe_impl.hpp l2hofe_impl.hpp \
diffop_impl.hpp


