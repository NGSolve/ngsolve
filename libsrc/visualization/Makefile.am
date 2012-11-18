noinst_HEADERS = meshdoc.hpp mvdraw.hpp vispar.hpp \
visual.hpp vssolution.hpp

include_HEADERS = soldata.hpp

AM_CPPFLAGS =  $(MPI_INCLUDES) -I$(top_srcdir)/libsrc/include  -DOPENGL -D$(TOGL_WINDOWINGSYSTEM) $(OCCFLAGS) $(TCL_INCLUDES)
METASOURCES = AUTO

if NGGUI
noinst_LTLIBRARIES = libvisual.la
endif


libvisual_la_SOURCES = meshdoc.cpp mvdraw.cpp \
	vsfieldlines.cpp vsmesh.cpp vssolution.cpp importsolution.cpp


