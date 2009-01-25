include_HEADERS = nglib.h

AM_CPPFLAGS = -I$(top_srcdir)/libsrc/include $(OCCFLAGS)  $(MPI_INCLUDES) 

lib_LTLIBRARIES = libnglib.la
libnglib_la_SOURCES = nglib.cpp

libnglib_la_LIBADD = \
	$(top_builddir)/libsrc/interface/libinterface.la \
	$(top_builddir)/libsrc/geom2d/libgeom2d.la \
	$(top_builddir)/libsrc/csg/libcsg.la \
	$(top_builddir)/libsrc/stlgeom/libstl.la \
	$(top_builddir)/libsrc/parallel/libparallel.la \
	$(top_builddir)/libsrc/meshing/libmesh.la \
	$(top_builddir)/libsrc/gprim/libgprim.la \
	$(top_builddir)/libsrc/opti/libopti.la \
	$(top_builddir)/libsrc/linalg/libla.la \
	$(top_builddir)/libsrc/general/libgeneral.la  $(MPI_LIBS)

libnglib_la_LDFLAGS = -avoid-version
#  -rdynamic


bin_PROGRAMS = ng_vol
ng_vol_SOURCES = ng_vol.cpp
ng_vol_LDADD = \
	libnglib.la 


#	$(top_builddir)/libsrc/occ/libocc.la 
#	$(top_builddir)/libsrc/geom2d/libgeom2d.la 
#	$(top_builddir)/libsrc/csg/libcsg.la 
#	$(top_builddir)/libsrc/interface/libinterface.la \
#	$(top_builddir)/libsrc/stlgeom/libstl.la \
#	$(top_builddir)/libsrc/meshing/libmesh.la 
#	$(top_builddir)/libsrc/gprim/libgprim.la 
#	$(top_builddir)/libsrc/opti/libopti.la 
#	$(top_builddir)/libsrc/linalg/libla.la \
#	$(top_builddir)/libsrc/general/libgeneral.la  


# ng_vol_LDFLAGS = -rdynamic



