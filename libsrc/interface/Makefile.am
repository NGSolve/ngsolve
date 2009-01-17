include_HEADERS = nginterface.h
noinst_HEADERS = writeuser.hpp

AM_CPPFLAGS = -I$(top_srcdir)/libsrc/include -I$(top_srcdir)/libsrc/interface $(TCL_INCLUDES) -DOPENGL
METASOURCES = AUTO
noinst_LIBRARIES = libinterface.a
libinterface_a_SOURCES = importsolution.cpp nginterface.cpp printdest.cpp \
	read_fnf_mesh.cpp readtetmesh.cpp readuser.cpp writeabaqus.cpp writediffpack.cpp \
	writedolfin.cpp writeelmer.cpp writefeap.cpp writefluent.cpp writegmsh.cpp writejcm.cpp \
	writepermas.cpp writetecplot.cpp writetet.cpp writetochnog.cpp writeuser.cpp \
	wuchemnitz.cpp


