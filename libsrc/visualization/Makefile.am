noinst_HEADERS = meshdoc.hpp mvdraw.hpp vispar.hpp \
visual.hpp vssolution.hpp

include_HEADERS = soldata.hpp

AM_CPPFLAGS =  $(MPI_INCLUDES) -I$(top_srcdir)/libsrc/include  -DOPENGL -D$(TOGL_WINDOWINGSYSTEM) $(OCCFLAGS) $(TCL_INCLUDES)
METASOURCES = AUTO

lib_LTLIBRARIES = libvisual.la


libvisual_la_SOURCES =

if NGGUI
libvisual_la_LIBADD =  $(LIBGLU) $(TCL_LIB_SPEC)
libvisual_la_SOURCES += meshdoc.cpp mvdraw.cpp \
	vsfieldlines.cpp vsmesh.cpp vssolution.cpp importsolution.cpp visualpkg.cpp
else
libvisual_la_SOURCES += visual_dummy.cpp
endif





