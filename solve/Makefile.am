AM_CPPFLAGS =  -I$(top_builddir) -I$(top_srcdir)/include $(PARDISO_INCLUDES) $(MUMPS_INCLUDES) $(SUPERLU_INCLUDES) $(TCL_INCLUDES) 

lib_LTLIBRARIES = libngsolve.la
bin_PROGRAMS =  ngs




libngsolve_la_SOURCES = bvp.cpp generalizedEVP.cpp \
numproc.cpp evp.cpp hyperbolic.cpp  \
numprocee.cpp pde.cpp pdeparser.cpp python_solve.cpp

if NGSGUI
  libngsolve_la_SOURCES += ngsolve.cpp shapetester.cpp
  AM_CPPFLAGS += -DTCL
endif

libngsolve_la_LIBADD = $(top_builddir)/comp/libngcomp.la	    \
$(top_builddir)/ngstd/libngstd.la 		\
$(top_builddir)/fem/libngfem.la 		\
$(top_builddir)/parallel/libparallel.la		\
$(top_builddir)/linalg/libngla.la		\
$(top_builddir)/basiclinalg/libngbla.la		\
$(top_builddir)/multigrid/libngmg.la		\
$(PARDISO_LIBS) $(LAPACK_LIBS) $(TCL_LIB_SPEC)

libngsolve_la_LDFLAGS = -avoid-version
# libngsolve_la_LDFLAGS += -lboost_python-py34 -lpython3.4m
# libngsolve_la_LDFLAGS += $(PYTHON_LIBS)
#   -rdynamic






if NGSPY
lib_LTLIBRARIES += libngspy.la
libngspy_la_SOURCES = ngspy.cpp bvp.cpp generalizedEVP.cpp \
numproc.cpp evp.cpp hyperbolic.cpp  \
numprocee.cpp pde.cpp pdeparser.cpp  python_solve.cpp

libngspy_la_LIBADD = $(top_builddir)/comp/libngcomp.la	    \
$(top_builddir)/ngstd/libngstd.la 		\
$(top_builddir)/fem/libngfem.la 		\
$(top_builddir)/parallel/libparallel.la		\
$(top_builddir)/linalg/libngla.la		\
$(top_builddir)/basiclinalg/libngbla.la		\
$(top_builddir)/multigrid/libngmg.la		\
$(PARDISO_LIBS) $(LAPACK_LIBS) $(TCL_LIB_SPEC)  \
-L$(libdir) -lnglib -linterface -lmesh -lstl -locc -lgeom2d -lcsg 
libngspy_la_LDFLAGS = -avoid-version
endif







ngs_SOURCES = ngs.cpp
ngs_LDADD = libngsolve.la \
	$(top_builddir)/comp/libngcomp.la \
	$(top_builddir)/multigrid/libngmg.la \
	$(top_builddir)/fem/libngfem.la \
	$(top_builddir)/parallel/libparallel.la \
	$(top_builddir)/linalg/libngla.la \
	$(top_builddir)/basiclinalg/libngbla.la		\
	$(top_builddir)/ngstd/libngstd.la \
	$(PARDISO_LIBS) $(LAPACK_LIBS) $(MUMPS_LIBS) \
	-L$(libdir) -lnglib -linterface -lmesh -lstl -locc -lgeom2d -lcsg 
#	-L/opt/OpenCASCADE/lib -lTKernel -lTKGeomBase -lTKMath -lTKG2d -lTKG3d -lTKXSBase -lTKOffset -lTKFillet -lTKShHealing -lTKMesh -lTKMeshVS -lTKTopAlgo -lTKGeomAlgo -lTKBool -lTKPrim -lTKBO -lTKIGES -lTKBRep -lTKSTEPBase -lTKSTEP -lTKSTL -lTKSTEPAttr -lTKSTEP209 -lTKXDESTEP -lTKXDEIGES -lTKXCAF -lTKLCAF -lFWOSPlugin




include_HEADERS = solve.hpp pde.hpp numproc.hpp 



