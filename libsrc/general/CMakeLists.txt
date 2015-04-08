add_definitions(-DNGINTERFACE_EXPORTS)
add_library(gen OBJECT
        array.cpp bitarray.cpp dynamicmem.cpp flags.cpp 
        hashtabl.cpp mystring.cpp ngexception.cpp optmem.cpp parthreads.cpp 
        profiler.cpp seti.cpp sort.cpp spbita2d.cpp symbolta.cpp table.cpp 
        mpi_interface.cpp gzstream.cpp
        )

set_target_properties( gen PROPERTIES POSITION_INDEPENDENT_CODE ON )

install( FILES ngexception.hpp DESTINATION include COMPONENT netgen_devel )

