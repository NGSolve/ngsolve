AM_CPPFLAGS = -I$(top_srcdir)/include $(NETGEN_INCLUDES)

lib_LTLIBRARIES = libngstd.la

libngstd_la_SOURCES = exception.cpp table.cpp bitarray.cpp flags.cpp   \
symboltable.cpp blockalloc.cpp evalfunc.cpp templates.cpp	       \
localheap.cpp stringops.cpp profiler.cpp array.hpp bessel.hpp	       \
exception.hpp memusage.hpp symboltable.hpp autodiffdiff.hpp	       \
bitarray.hpp flags.hpp ngstd.hpp table.hpp autodiff.hpp blockalloc.hpp \
hashtable.hpp profiler.hpp templates.hpp autoptr.hpp evalfunc.hpp      \
localheap.hpp stringops.hpp


libngstd_la_LDFLAGS = -avoid-version