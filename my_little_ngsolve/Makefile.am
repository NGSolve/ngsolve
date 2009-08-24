AM_CPPFLAGS = -I$(top_srcdir)/include $(NETGEN_INCLUDES) $(LAPACK_INCLUDES) $(TCL_INCLUDES)

lib_LTLIBRARIES = libmyngsolve.la

libmyngsolve_la_SOURCES = myElement.cpp myFESpace.cpp myIntegrator.cpp \
myElement.hpp myFESpace.hpp myIntegrator.hpp myHOElement.hpp	       \
myHOElement.cpp myHOFESpace.hpp myHOFESpace.cpp all_in_one.cpp         \
demo_instat.cpp demo_stokes.cpp demo_coupling.cpp demo_nonlinear.cpp   

dist_noinst_DATA = simple.pde all_in_one.pde square.in2d square.vol \
stokes.pde instat.pde nonlinear.pde coupling.pde



