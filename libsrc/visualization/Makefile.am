noinst_HEADERS = meshdoc.hpp mvdraw.hpp vispar.hpp \
visual.hpp vssolution.hpp

include_HEADERS = soldata.hpp

AM_CPPFLAGS = -I$(top_srcdir)/libsrc/include  -DOPENGL  $(OCCFLAGS)    $(TCL_INCLUDES)
METASOURCES = AUTO
noinst_LIBRARIES = libvisual.a
libvisual_a_SOURCES = meshdoc.cpp mvdraw.cpp stlmeshing.cpp vscsg.cpp \
	vsfieldlines.cpp vsmesh.cpp vsocc.cpp vssolution.cpp
AM_CXXFLAGS = -DOPENGL
