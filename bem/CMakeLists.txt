include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})

add_library( ngsbem ${NGS_LIB_TYPE}
        ${ngsbem_object_libs}
        mptools.cpp potentialtools.cpp python_bem.cpp ngbem.cpp
        intrules_SauterSchwab.cpp analytic_integrals.cpp
        )

target_include_directories(ngsbem PRIVATE ${NETGEN_PYTHON_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}/../ngstd ${CMAKE_CURRENT_SOURCE_DIR}/../linalg)

target_link_libraries (ngsbem PUBLIC ngcomp ngfem ngla ngbla ngstd ${MPI_CXX_LIBRARIES} PRIVATE "$<BUILD_INTERFACE:netgen_python>" ${HYPRE_LIBRARIES})
target_link_libraries(ngsbem ${LAPACK_CMAKE_LINK_INTERFACE} "$<BUILD_INTERFACE:ngs_lapack>")
install( TARGETS ngsbem ${ngs_install_dir} )

install( FILES
        mptools.hpp potentialtools.hpp mp_coefficient.hpp
        intrules_SauterSchwab.hpp analytic_integrals.hpp
        DESTINATION ${NGSOLVE_INSTALL_DIR_INCLUDE}
        COMPONENT ngsolve_devel
       )


