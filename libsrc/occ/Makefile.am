
AM_CPPFLAGS = -I$(top_srcdir)/libsrc/include  $(OCCFLAGS) 

# $(OCC_INC_FLAG)

METASOURCES = AUTO

noinst_LIBRARIES = libocc.a

libocc_a_SOURCES = Partition_Inter2d.cxx Partition_Inter3d.cxx \
	Partition_Loop.cxx Partition_Loop2d.cxx Partition_Loop3d.cxx Partition_Spliter.cxx \
	occconstruction.cpp occgenmesh.cpp occgeom.cpp occmeshsurf.cpp
