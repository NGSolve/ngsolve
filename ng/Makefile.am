noinst_HEADERS = demoview.hpp

AM_CPPFLAGS = -I$(top_srcdir)/libsrc/include -I$(top_srcdir)/libsrc/interface -DOPENGL $(OCCFLAGS)  $(TCL_INCLUDES) 

bin_PROGRAMS = netgen
netgen_SOURCES = demoview.cpp ngappinit.cpp ngpkg.cpp onetcl.cpp


netgen_LDADD = $(top_builddir)/libsrc/csg/libcsg.a \
	$(top_builddir)/libsrc/geom2d/libgeom2d.a $(top_builddir)/libsrc/interface/libinterface.a \
	$(top_builddir)/libsrc/stlgeom/libstl.a $(top_builddir)/libsrc/occ/libocc.a \
	$(top_builddir)/libsrc/visualization/libvisual.a $(top_builddir)/libsrc/meshing/libmesh.a \
	$(top_builddir)/libsrc/gprim/libgprim.a $(top_builddir)/libsrc/opti/libopti.a \
	$(top_builddir)/libsrc/linalg/libla.a $(top_builddir)/libsrc/general/libgeneral.a \
	$(OCCLIBS) $(TOGLLIBDIR) -lTogl1.7 -lGLU $(TK_LIB_SPEC) $(TCL_LIB_SPEC)




dist_bin_SCRIPTS = dialog.tcl menustat.tcl ngicon.tcl ng.tcl \
ngvisual.tcl sockets.tcl drawing.tcl nghelp.tcl ngshell.tcl  \
ngtesting.tcl parameters.tcl variables.tcl


netgen_LDFLAGS = -rdynamic
