noinst_HEADERS = meshstlsurface.hpp stlgeom.hpp stlline.hpp \
stltool.hpp stltopology.hpp vsstl.hpp

AM_CPPFLAGS = -I$(top_srcdir)/libsrc/include $(TCL_INCLUDES)
METASOURCES = AUTO

lib_LTLIBRARIES = libstl.la 


libstl_la_SOURCES = meshstlsurface.cpp stlgeom.cpp stlgeomchart.cpp \
	stlgeommesh.cpp stlline.cpp stltool.cpp stltopology.cpp


libstl_la_LIBADD =  $(top_builddir)/libsrc/meshing/libmesh.la




if NGGUI
lib_LTLIBRARIES += libstlvis.la 
libstlvis_la_SOURCES = stlpkg.cpp 
libstlvis_la_LIBADD = libstl.la  $(top_builddir)/libsrc/visualization/libvisual.la
libstl_la_SOURCES += vsstl.cpp
libstl_la_LIBADD += $(top_builddir)/libsrc/visualization/libvisual.la
endif
