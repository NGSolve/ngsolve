noinst_HEADERS = densemat.hpp linalg.hpp polynomial.hpp vector.hpp opti.hpp
AM_CPPFLAGS = -I$(top_srcdir)/libsrc/include
METASOURCES = AUTO
noinst_LTLIBRARIES = libla.la
libla_la_SOURCES = densemat.cpp polynomial.cpp bfgs.cpp linopt.cpp linsearch.cpp

#  vector.cpp

libla_la_LDFLAGS = -rdynamic
