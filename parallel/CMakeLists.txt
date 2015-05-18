include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})

add_library(parallel OBJECT
        parallelvvector.cpp parallel_matrices.cpp 
        )

set_target_properties( parallel PROPERTIES POSITION_INDEPENDENT_CODE ON )

install( FILES
        parallelngs.hpp parallelvector.hpp parallel_matrices.hpp dump.hpp
        DESTINATION include
        COMPONENT ngsolve_devel
       )

