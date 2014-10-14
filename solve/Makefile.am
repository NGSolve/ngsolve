AM_CPPFLAGS =  -I$(top_builddir) -I$(top_srcdir)/include $(PARDISO_INCLUDES) $(MUMPS_INCLUDES) $(SUPERLU_INCLUDES) $(TCL_INCLUDES) 

lib_LTLIBRARIES = libsolve.la
bin_PROGRAMS =  ngs


libsolve_la_SOURCES = bvp.cpp generalizedEVP.cpp numproc.cpp evp.cpp \
hyperbolic.cpp numprocee.cpp python_solve.cpp python_mpi.cpp


libsolve_la_LIBADD = $(top_builddir)/comp/libngcomp.la		    \
$(top_builddir)/ngstd/libngstd.la $(top_builddir)/fem/libngfem.la   \
$(top_builddir)/parallel/libparallel.la				    \
$(top_builddir)/linalg/libngla.la				    \
$(top_builddir)/basiclinalg/libngbla.la				    \
$(top_builddir)/multigrid/libngmg.la $(PARDISO_LIBS) $(LAPACK_LIBS) 

libsolve_la_LDFLAGS = " -avoid-version -rdynamic"


if NGSGUI
lib_LTLIBRARIES += libngsolve.la	
  libngsolve_la_SOURCES = ngsolve.cpp shapetester.cpp
  AM_CPPFLAGS += -DTCL
  libngsolve_la_LDFLAGS = " -avoid-version -rdynamic"

  libngsolve_la_LIBADD = libsolve.la			   \
	$(top_builddir)/comp/libngcomp.la		   \
	$(top_builddir)/ngstd/libngstd.la		   \
	$(top_builddir)/fem/libngfem.la			   \
	$(top_builddir)/parallel/libparallel.la		   \
	$(top_builddir)/linalg/libngla.la		   \
	$(top_builddir)/basiclinalg/libngbla.la		   \
	$(top_builddir)/multigrid/libngmg.la		   \
	-L$(libdir) -linterface				   \
	$(PARDISO_LIBS) $(LAPACK_LIBS) $(TCL_LIB_SPEC)	  
endif



ngs_SOURCES = ngs.cpp
ngs_LDADD = libsolve.la \
	$(top_builddir)/comp/libngcomp.la \
	$(top_builddir)/multigrid/libngmg.la \
	$(top_builddir)/fem/libngfem.la \
	$(top_builddir)/parallel/libparallel.la \
	$(top_builddir)/linalg/libngla.la \
	$(top_builddir)/basiclinalg/libngbla.la		\
	$(top_builddir)/ngstd/libngstd.la \
	$(PARDISO_LIBS) $(LAPACK_LIBS) $(MUMPS_LIBS) \
	-L$(libdir) -linterface


include_HEADERS = solve.hpp



