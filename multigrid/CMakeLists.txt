include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})

install( FILES
        mgpre.hpp  multigrid.hpp  prolongation.hpp  smoother.hpp
        DESTINATION ${NGSOLVE_INSTALL_DIR_INCLUDE}
        COMPONENT ngsolve_devel
       )

