bin_PROGRAMS = demo_std demo_bla demo_fem 
# demo_comp

AM_CPPFLAGS = -I$(top_srcdir)/include -DLAPACK

demo_std_SOURCES = demo_std.cpp
demo_std_LDADD = $(top_builddir)/ngstd/libngstd.la $(MPI_LIBS)


demo_bla_SOURCES = demo_bla.cpp
demo_bla_LDADD = $(top_builddir)/basiclinalg/libngbla.la \
	$(top_builddir)/ngstd/libngstd.la $(LAPACK_LIBS) $(MPI_LIBS)


demo_fem_SOURCES = demo_fem.cpp
demo_fem_LDADD = $(top_builddir)/fem/libngfem.la \
	$(top_builddir)/basiclinalg/libngbla.la \
	$(top_builddir)/ngstd/libngstd.la $(LAPACK_LIBS) $(MPI_LIBS)


# demo_comp_SOURCES = demo_comp.cpp
#demo_comp_LDADD = \
#	$(top_builddir)/comp/libngcomp.la \
#	$(top_builddir)/fem/libngfem.la \
#	$(top_builddir)/linalg/libngla.la \
#	$(top_builddir)/ngstd/libngstd.la \
#	$(top_builddir)/solve/libngsolve.la \
#	$(top_builddir)/parallel/libparallel.la \
#	-L/opt/netgen/lib -lnglib -linterface $(LAPACK_LIBS) 
# demo_comp_CPPFLAGS = $(NETGEN_INCLUDES)  -I$(top_srcdir)/include




