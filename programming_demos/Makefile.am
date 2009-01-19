bin_PROGRAMS = demo_std demo_bla demo_fem

AM_CPPFLAGS = -I$(top_srcdir)/include

demo_std_SOURCES = demo_std.cpp
demo_std_LDADD = $(top_builddir)/ngstd/libngstd.la


demo_bla_SOURCES = demo_bla.cpp
demo_bla_LDADD = $(top_builddir)/basiclinalg/libngbla.la \
	$(top_builddir)/ngstd/libngstd.la $(LAPACK_LIB)


demo_fem_SOURCES = demo_fem.cpp
demo_fem_LDADD = $(top_builddir)/fem/libngfem.la \
	$(top_builddir)/basiclinalg/libngbla.la \
	$(top_builddir)/ngstd/libngstd.la $(LAPACK_LIB)
