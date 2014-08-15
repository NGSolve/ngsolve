AM_CPPFLAGS = -I$(top_builddir) -I$(top_srcdir)/include
#  -I/usr/include/python$(PYTHON_VERSION)

lib_LTLIBRARIES = libngstd.la

include_HEADERS = array.hpp autodiffdiff.hpp autodiff.hpp autoptr.hpp \
bessel.hpp bitarray.hpp blockalloc.hpp evalfunc.hpp exception.hpp     \
flags.hpp hashtable.hpp localheap.hpp memusage.hpp ngstd.hpp	      \
profiler.hpp stringops.hpp symboltable.hpp table.hpp templates.hpp    \
parthreads.hpp statushandler.hpp ngsstream.hpp mpiwrapper.hpp	      \
polorder.hpp archive.hpp archive_base.hpp sockets.hpp cuda_ngstd.hpp  \
mycomplex.hpp python_ngstd.hpp


libngstd_la_SOURCES = exception.cpp table.cpp bitarray.cpp flags.cpp \
symboltable.cpp blockalloc.cpp evalfunc.cpp templates.cpp	         \
localheap.cpp stringops.cpp profiler.cpp archive.cpp sockets.cpp     \
cuda_ngstd.cpp python_ngstd.cpp

# why are the headers in sources ???? JS
# array.hpp bessel.hpp exception.hpp memusage.hpp symboltable.hpp	     \
# autodiffdiff.hpp bitarray.hpp flags.hpp ngstd.hpp table.hpp	     \
# autodiff.hpp blockalloc.hpp hashtable.hpp profiler.hpp templates.hpp \
# autoptr.hpp evalfunc.hpp localheap.hpp stringops.hpp mpiwrapper.hpp  \
# polorder.hpp 


libngstd_la_LDFLAGS = -avoid-version

dist_bin_SCRIPTS = init.py expr.py
