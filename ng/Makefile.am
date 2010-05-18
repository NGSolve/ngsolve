include_HEADERS = 

AM_CPPFLAGS = -I$(top_srcdir)/libsrc/include -I$(top_srcdir)/libsrc/interface -DOPENGL  -D$(TOGL_WINDOWINGSYSTEM) $(OCCFLAGS)  $(TCL_INCLUDES) $(MPI_INCLUDES) $(FFMPEG_INCLUDES) $(JPEGLIB_INCLUDES)

bin_PROGRAMS = netgen
netgen_SOURCES = demoview.cpp ngappinit.cpp ngpkg.cpp onetcl.cpp nginterface.cpp nginterface_v2.cpp parallelfunc.cpp parallelinterface.cpp demoview.hpp parallelfunc.hpp togl_1_7.h


netgen_LDADD = $(top_builddir)/libsrc/visualization/libvisual.a \
	$(top_builddir)/libsrc/csg/libcsg.la \
	$(top_builddir)/libsrc/geom2d/libgeom2d.la \
	$(top_builddir)/libsrc/interface/libinterface.la \
	$(top_builddir)/libsrc/stlgeom/libstl.la \
	$(top_builddir)/libsrc/occ/libocc.la \
	$(top_builddir)/libsrc/meshing/libmesh.la \
	$(top_builddir)/libsrc/gprim/libgprim.la \
	$(top_builddir)/libsrc/linalg/libla.la \
	$(top_builddir)/libsrc/general/libgen.la \
	$(OCCLIBS) -L$(TK_BIN_DIR)/Togl1.7 $(TOGLLIBDIR) -lTogl1.7 $(LIBGLU) $(TK_LIB_SPEC) $(TCL_LIB_SPEC) $(MPI_LIBS) $(FFMPEG_LIBS) $(JPEGLIB_LIBS) $(PKG_LIBS) 

# add for static linkage of ngsolve:
#   /opt/netgen/lib/libngsolve.a /opt/netgen/lib/libngcomp.a /opt/netgen/lib/libngcomp.a  /opt/netgen/lib/libngfemng.a   /opt/netgen/lib/libngmg.a  /opt/netgen/lib/libngla.a  /opt/netgen/lib/libngbla.a  /opt/netgen/lib/libngstd.a -L/opt/intel/mkl/10.2.1.017/lib/em64t /opt/intel/mkl/10.2.1.017/lib/em64t/libmkl_solver_lp64.a  -lmkl_intel_lp64  -lmkl_gnu_thread -lmkl_core
#
#

dist_bin_SCRIPTS = dialog.tcl menustat.tcl ngicon.tcl ng.tcl \
ngvisual.tcl sockets.tcl drawing.tcl nghelp.tcl ngshell.tcl  \
ngtesting.tcl parameters.tcl variables.tcl netgen.ocf

netgen_LDFLAGS = -rdynamic  -static


