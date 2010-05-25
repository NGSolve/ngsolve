include_HEADERS = nglib.h

dist_pkgdata_DATA = cube.surf

AM_CPPFLAGS = -I$(top_srcdir)/libsrc/include  $(MPI_INCLUDES) $(OCCFLAGS)

lib_LTLIBRARIES = libnglib.la
libnglib_la_SOURCES = nglib.cpp

libnglib_la_LIBADD = \
	$(top_builddir)/libsrc/interface/libinterface.la \
	$(top_builddir)/libsrc/geom2d/libgeom2d.la \
	$(top_builddir)/libsrc/csg/libcsg.la \
	$(top_builddir)/libsrc/stlgeom/libstl.la \
	$(top_builddir)/libsrc/occ/libocc.la \
	$(top_builddir)/libsrc/meshing/libmesh.la \
	$(top_builddir)/libsrc/gprim/libgprim.la \
	$(top_builddir)/libsrc/linalg/libla.la \
	$(top_builddir)/libsrc/general/libgen.la \
	$(OCCLIBS) $(MPI_LIBS)

libnglib_la_LDFLAGS = -avoid-version
#  -rdynamic


bin_PROGRAMS = ng_vol ng_stl
#  ng_occ

ng_vol_SOURCES = ng_vol.cpp
ng_vol_LDADD = \
	libnglib.la 
ng_stl_SOURCES = ng_stl.cpp
ng_stl_LDADD = \
	libnglib.la 


# ng_occ_SOURCES = ng_occ.cpp
# ng_occ_LDADD = \
#	libnglib.la 


