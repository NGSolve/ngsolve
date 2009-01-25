include_HEADERS = 

AM_CPPFLAGS = -I$(top_srcdir)/libsrc/include -I$(top_srcdir)/libsrc/interface -DOPENGL $(OCCFLAGS)  $(TCL_INCLUDES) $(MPI_INCLUDES) 

bin_PROGRAMS = netgen
netgen_SOURCES = demoview.cpp ngappinit.cpp ngpkg.cpp onetcl.cpp nginterface.cpp parallelfunc.cpp parallelinterface.cpp demoview.hpp parallelfunc.hpp togl_1_7.h


netgen_LDADD = $(top_builddir)/libsrc/visualization/libvisual.a \
	$(top_builddir)/libsrc/csg/libcsg.la \
	$(top_builddir)/libsrc/geom2d/libgeom2d.la \
	$(top_builddir)/libsrc/interface/libinterface.la \
	$(top_builddir)/libsrc/stlgeom/libstl.la \
	$(top_builddir)/libsrc/occ/libocc.la \
	$(top_builddir)/libsrc/meshing/libmesh.la \
	$(top_builddir)/libsrc/gprim/libgprim.la \
	$(top_builddir)/libsrc/linalg/libla.la \
	$(top_builddir)/libsrc/general/libgeneral.la \
	$(OCCLIBS) -L$(TK_BIN_DIR)/Togl1.7 $(TOGLLIBDIR) -lTogl1.7 -lGLU $(TK_LIB_SPEC) $(TCL_LIB_SPEC) $(MPI_LIBS)




dist_bin_SCRIPTS = dialog.tcl menustat.tcl ngicon.tcl ng.tcl \
ngvisual.tcl sockets.tcl drawing.tcl nghelp.tcl ngshell.tcl  \
ngtesting.tcl parameters.tcl variables.tcl


netgen_LDFLAGS = -rdynamic
