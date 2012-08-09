bin_PROGRAMS = demo_std demo_bla demo_fem demo_comp demo_comp1d demo_solve

AM_CPPFLAGS = -I$(top_srcdir)/include -DLAPACK  $(LAPACK_INCLUDES) $(MPI_INCLUDES)

demo_std_SOURCES = demo_std.cpp
demo_std_LDADD = $(top_builddir)/ngstd/libngstd.la $(MPI_LIBS)


demo_bla_SOURCES = demo_bla.cpp
demo_bla_LDADD = $(top_builddir)/basiclinalg/libngbla.la \
	$(top_builddir)/ngstd/libngstd.la $(LAPACK_LIBS) $(MPI_LIBS)


demo_fem_SOURCES = demo_fem.cpp
demo_fem_LDADD = $(top_builddir)/fem/libngfem.la \
	$(top_builddir)/basiclinalg/libngbla.la \
	$(top_builddir)/ngstd/libngstd.la $(LAPACK_LIBS) $(MPI_LIBS)


demo_comp_SOURCES = demo_comp.cpp
demo_comp_LDADD = \
	$(top_builddir)/comp/libngcomp.la \
	$(top_builddir)/solve/libngsolve.la \
	$(top_builddir)/fem/libngfem.la \
	$(top_builddir)/linalg/libngla.la \
	$(top_builddir)/ngstd/libngstd.la \
	$(top_builddir)/parallel/libparallel.la \
	$(LAPACK_LIBS) \
        -L$(libdir) -lnglib -linterface 
demo_comp_CPPFLAGS = -I$(top_srcdir)/include $(NETGEN_INCLUDES) $(LAPACK_INCLUDES) $(MPI_INCLUDES)


demo_comp1d_SOURCES = demo_comp1d.cpp netgen_dummy.cpp
demo_comp1d_LDADD = \
	$(top_builddir)/comp/libngcomp.la \
	$(top_builddir)/solve/libngsolve.la \
	$(top_builddir)/fem/libngfem.la \
	$(top_builddir)/linalg/libngla.la \
	$(top_builddir)/ngstd/libngstd.la \
	$(top_builddir)/parallel/libparallel.la \
	$(LAPACK_LIBS) 
#        -L$(libdir) -lnglib -linterface 
demo_comp1d_CPPFLAGS = -I$(top_srcdir)/include $(NETGEN_INCLUDES) $(LAPACK_INCLUDES) $(MPI_INCLUDES)




demo_solve_SOURCES = demo_solve.cpp
demo_solve_LDADD = \
	$(top_builddir)/solve/libngsolve.la \
	$(top_builddir)/comp/libngcomp.la \
	$(top_builddir)/fem/libngfem.la \
	$(top_builddir)/parallel/libparallel.la \
	$(top_builddir)/linalg/libngla.la \
	$(top_builddir)/ngstd/libngstd.la \
	$(LAPACK_LIBS) \
        -L$(libdir) -lnglib -linterface 
demo_solve_CPPFLAGS = -I$(top_srcdir)/include $(NETGEN_INCLUDES) $(LAPACK_INCLUDES) $(MPI_INCLUDES)


