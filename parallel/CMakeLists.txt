include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})

add_library(parallel OBJECT
        parallelvvector.cpp parallel_matrices.cpp 
        )

set_target_properties( parallel PROPERTIES POSITION_INDEPENDENT_CODE ON )
target_compile_definitions(parallel PUBLIC ${NGSOLVE_COMPILE_DEFINITIONS})
target_compile_options(parallel PUBLIC ${NGSOLVE_COMPILE_OPTIONS})
target_include_directories(parallel PUBLIC ${NETGEN_INCLUDE_DIRS})

install( FILES
        parallelngs.hpp parallelvector.hpp parallel_matrices.hpp dump.hpp
        DESTINATION ${NGSOLVE_INSTALL_DIR_INCLUDE}
        COMPONENT ngsolve_devel
       )

