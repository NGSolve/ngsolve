noinst_HEADERS = writeuser.hpp

AM_CPPFLAGS = -I$(top_srcdir)/libsrc/include -I$(top_srcdir)/libsrc/interface  $(MPI_INCLUDES) $(TCL_INCLUDES) -DOPENGL
METASOURCES = AUTO
lib_LTLIBRARIES = libinterface.la
libinterface_la_SOURCES = nginterface.cpp nginterface_v2.cpp \
	read_fnf_mesh.cpp readtetmesh.cpp readuser.cpp writeabaqus.cpp writediffpack.cpp \
	writedolfin.cpp writeelmer.cpp writefeap.cpp writefluent.cpp writegmsh.cpp writejcm.cpp \
	writepermas.cpp writetecplot.cpp writetet.cpp writetochnog.cpp writeuser.cpp \
	wuchemnitz.cpp writegmsh2.cpp writeOpenFOAM15x.cpp 


libinterface_la_LIBADD = $(top_builddir)/libsrc/meshing/libmesh.la \
	$(top_builddir)/libsrc/visualization/libvisual.la \
	$(top_builddir)/libsrc/csg/libcsg.la \
	$(top_builddir)/libsrc/geom2d/libgeom2d.la

# libinterface_la_LDFLAGS = -rdynamic
