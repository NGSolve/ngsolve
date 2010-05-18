noinst_HEADERS = array.hpp myadt.hpp optmem.hpp sort.hpp table.hpp autodiff.hpp flags.hpp mystring.hpp spbita2d.hpp template.hpp autoptr.hpp hashtabl.hpp netgenout.hpp profiler.hpp stack.hpp bitarray.hpp seti.hpp symbolta.hpp dynamicmem.hpp  parthreads.hpp mpi_interface.hpp

#  moveablemem.hpp

include_HEADERS = ngexception.hpp

AM_CPPFLAGS =  $(MPI_INCLUDES) -I$(top_srcdir)/libsrc/include
METASOURCES = AUTO
noinst_LTLIBRARIES = libgen.la
libgen_la_SOURCES = array.cpp bitarray.cpp dynamicmem.cpp flags.cpp \
	hashtabl.cpp mystring.cpp ngexception.cpp optmem.cpp parthreads.cpp \
	profiler.cpp seti.cpp sort.cpp spbita2d.cpp symbolta.cpp table.cpp

#  moveablemem.cpp
