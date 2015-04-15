include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_library(ngmg OBJECT
        mgpre.cpp  prolongation.cpp  smoother.cpp
        )

set_target_properties( ngmg PROPERTIES POSITION_INDEPENDENT_CODE ON )

install( FILES
        mgpre.hpp  multigrid.hpp  prolongation.hpp  smoother.hpp
        DESTINATION include
        COMPONENT ngsolve_devel
       )

