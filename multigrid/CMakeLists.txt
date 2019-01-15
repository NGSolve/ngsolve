include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})

add_library(ngmg INTERFACE)
set(sdir ${CMAKE_CURRENT_SOURCE_DIR})
target_sources(ngmg INTERFACE ${sdir}/mgpre.cpp ${sdir}/prolongation.cpp ${sdir}/smoother.cpp)
target_link_libraries(ngmg INTERFACE ngcore)

install( FILES
        mgpre.hpp  multigrid.hpp  prolongation.hpp  smoother.hpp
        DESTINATION ${NGSOLVE_INSTALL_DIR_INCLUDE}
        COMPONENT ngsolve_devel
       )

