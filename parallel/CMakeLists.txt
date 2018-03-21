include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})

install( FILES
        parallelngs.hpp parallelvector.hpp parallel_matrices.hpp dump.hpp
        DESTINATION ${NGSOLVE_INSTALL_DIR_INCLUDE}
        COMPONENT ngsolve_devel
       )

