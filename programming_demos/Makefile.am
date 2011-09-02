bin_PROGRAMS = demo_std demo_bla demo_fem demo_comp
# demo_comp

AM_CPPFLAGS = -I$(top_srcdir)/include

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
demo_comp_LDADD =  $(top_builddir)/solve/libngsolve.la \
	 -L/opt/netgen/lib -lnglib $(LAPACK_LIBS) 
#  $(TCL_LIB_SPEC)
#
# $(top_builddir)/comp/libngcomp.la \
#	$(top_builddir)/comp/libngcomp.la $(LAPACK_LIBS) 

demo_comp_CPPFLAGS = -DNETGEN_ELTRANS -I/opt/netgen/include -I/home/joachim/netgen-mesher/netgen/libsrc/include -I/home/joachim/netgen-mesher/netgen/libsrc/interface  -I$(top_srcdir)/include

