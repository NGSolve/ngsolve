include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})

add_library(ngmg OBJECT
        mgpre.cpp  prolongation.cpp  smoother.cpp
        )

set_target_properties( ngmg PROPERTIES POSITION_INDEPENDENT_CODE ON )
target_link_libraries( ngmg PUBLIC ngs_interface PRIVATE ngs_private_interface)

install( FILES
        mgpre.hpp  multigrid.hpp  prolongation.hpp  smoother.hpp
        DESTINATION ${NGSOLVE_INSTALL_DIR_INCLUDE}
        COMPONENT ngsolve_devel
       )

