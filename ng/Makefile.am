
include_HEADERS = 

AM_CPPFLAGS = -I$(top_srcdir)/libsrc/include -DOPENGL  -D$(TOGL_WINDOWINGSYSTEM)  $(TCL_INCLUDES) $(MPI_INCLUDES) $(FFMPEG_INCLUDES) $(JPEGLIB_INCLUDES) 

# SUBDIRS = Togl-1.7

if NGGUI
bin_PROGRAMS = netgen
endif

netgen_SOURCES = demoview.cpp ngappinit.cpp onetcl.cpp parallelfunc.cpp ngpkg.cpp demoview.hpp parallelfunc.hpp togl_1_7.h
#  nginterface.cpp nginterface_v2.cpp


if NGMKL
AM_CPPFLAGS += -DLINKMKL
MKL_LIBS = -L$(MKLROOT)/lib/intel64  -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core -lmkl_def
endif
# -lpardiso500-GNU472-X86-64 -lgfortran

netgen_LDADD = $(top_builddir)/libsrc/visualization/libvisual.la \
	$(top_builddir)/libsrc/csg/libcsgvis.la \
	$(top_builddir)/libsrc/csg/libcsg.la \
	$(top_builddir)/libsrc/stlgeom/libstlvis.la \
	$(top_builddir)/libsrc/stlgeom/libstl.la \
	$(top_builddir)/libsrc/geom2d/libgeom2d.la \
	$(top_builddir)/libsrc/geom2d/libgeom2dvis.la \
	$(top_builddir)/libsrc/interface/libinterface.la \
	$(top_builddir)/libsrc/meshing/libmesh.la \
	-L$(TK_BIN_DIR)/../Togl1.7 -lTogl1.7 $(TOGLLIBDIR)  $(LIBGLU) $(TK_LIB_SPEC) $(TCL_LIB_SPEC) $(MPI_LIBS) $(FFMPEG_LIBS) $(JPEGLIB_LIBS) $(PKG_LIBS) $(MKL_LIBS) 
# -lTogl1.7

#
#	$(top_builddir)/libsrc/occ/liboccvis.la 
#	$(top_builddir)/libsrc/occ/libocc.la 
#	$(top_builddir)/libsrc/stlgeom/libstlvis.la 
#	$(top_builddir)/libsrc/stlgeom/libstl.la 
#	$(top_builddir)/libsrc/geom2d/libgeom2d.la 
#	$(top_builddir)/libsrc/geom2d/libgeom2dvis.la 


# add for static linkage of ngsolve:
#   /opt/netgen/lib/libngsolve.a /opt/netgen/lib/libngcomp.a /opt/netgen/lib/libngcomp.a  /opt/netgen/lib/libngfemng.a   /opt/netgen/lib/libngmg.a  /opt/netgen/lib/libngla.a  /opt/netgen/lib/libngbla.a  /opt/netgen/lib/libngstd.a -L/opt/intel/mkl/10.2.1.017/lib/em64t /opt/intel/mkl/10.2.1.017/lib/em64t/libmkl_solver_lp64.a  -lmkl_intel_lp64  -lmkl_gnu_thread -lmkl_core 
#
#



lib_LTLIBRARIES = libngpy.la	
libngpy_la_SOURCES = netgenpy.cpp
libngpy_la_LDFLAGS = -avoid-version -module

libngpy_la_LIBADD =  $(top_builddir)/libsrc/visualization/libvisual.la \
	$(top_builddir)/libsrc/csg/libcsg.la \
	$(top_builddir)/libsrc/stlgeom/libstl.la \
	$(top_builddir)/libsrc/geom2d/libgeom2d.la \
	$(top_builddir)/libsrc/interface/libinterface.la \
	$(top_builddir)/libsrc/meshing/libmesh.la



dist_bin_SCRIPTS = dialog.tcl menustat.tcl ngicon.tcl ng.tcl	  \
ngvisual.tcl sockets.tcl drawing.tcl nghelp.tcl ngshell.tcl	  \
ngtesting.tcl parameters.tcl variables.tcl csgeom.tcl stlgeom.tcl \
occgeom.tcl acisgeom.tcl netgen.ocf


netgen_LDFLAGS = -export-dynamic

# if NGPYTHON
# AM_CPPFLAGS += -I/usr/include/python3.4 -DNGPYTHON
# netgen_LDFLAGS += -lboost_python-py34 -lpython3.4m
# netgen_LDFLAGS += -lpython3.4m
# netgen_LDFLAGS += -lboost_python-py32 -lpython3.2mu	
# endif


# netgen_LDFLAGS = -rdynamic 
# -static

# netgen_LDFLAGS += -rdynamic



