noinst_HEADERS = meshstlsurface.hpp stlgeom.hpp stlline.hpp \
stltool.hpp stltopology.hpp vsstl.hpp

AM_CPPFLAGS = -I$(top_srcdir)/libsrc/include $(TCL_INCLUDES)
METASOURCES = AUTO

lib_LTLIBRARIES = libstl.la libstlvis.la

libstl_la_SOURCES = meshstlsurface.cpp stlgeom.cpp stlgeomchart.cpp \
	stlgeommesh.cpp stlline.cpp stltool.cpp stltopology.cpp


libstlvis_la_SOURCES = stlpkg.cpp vsstl.cpp
libstlvis_la_LIBADD = libstl.la $(top_builddir)/libsrc/linalg/libla.la 

