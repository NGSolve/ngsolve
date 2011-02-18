noinst_HEADERS = geom2dmesh.hpp  splinegeometry.hpp  spline.hpp geometry2d.hpp vsgeom2d.hpp

AM_CPPFLAGS = -I$(top_srcdir)/libsrc/include  $(TCL_INCLUDES)

METASOURCES = AUTO
noinst_LTLIBRARIES = libgeom2d.la libgeom2dvis.la

libgeom2d_la_SOURCES = genmesh2d.cpp geom2dmesh.cpp spline.cpp \
	splinegeometry.cpp


libgeom2dvis_la_SOURCES = geom2dpkg.cpp vsgeom2d.cpp
libgeom2dvis_la_LIBADD = libgeom2d.la


