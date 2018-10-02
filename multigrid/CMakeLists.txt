include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})

add_library(ngmg INTERFACE)
target_sources(ngmg INTERFACE
  ${CMAKE_CURRENT_SOURCE_DIR}/mgpre.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/prolongation.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/smoother.cpp
)

install( FILES
        mgpre.hpp  multigrid.hpp  prolongation.hpp  smoother.hpp
        DESTINATION ${NGSOLVE_INSTALL_DIR_INCLUDE}
        COMPONENT ngsolve_devel
       )

