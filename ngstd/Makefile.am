AM_CPPFLAGS = -I$(top_srcdir)/include

lib_LTLIBRARIES = libngstd.la

include_HEADERS = array.hpp autodiffdiff.hpp autodiff.hpp autoptr.hpp \
bessel.hpp bitarray.hpp blockalloc.hpp evalfunc.hpp exception.hpp     \
flags.hpp hashtable.hpp localheap.hpp memusage.hpp ngstd.hpp	      \
profiler.hpp stringops.hpp symboltable.hpp table.hpp templates.hpp    \
parthreads.hpp statushandler.hpp ngsstream.hpp mpiwrapper.hpp



libngstd_la_SOURCES = exception.cpp table.cpp bitarray.cpp flags.cpp   \
symboltable.cpp blockalloc.cpp evalfunc.cpp templates.cpp	       \
localheap.cpp stringops.cpp profiler.cpp array.hpp bessel.hpp	       \
exception.hpp memusage.hpp symboltable.hpp autodiffdiff.hpp	       \
bitarray.hpp flags.hpp ngstd.hpp table.hpp autodiff.hpp blockalloc.hpp \
hashtable.hpp profiler.hpp templates.hpp autoptr.hpp evalfunc.hpp      \
localheap.hpp stringops.hpp mpiwrapper.hpp


libngstd_la_LDFLAGS = -avoid-version