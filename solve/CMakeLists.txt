include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})

add_library( solve ${NGS_LIB_TYPE}
        python_solve.cpp
        )
#  python_mpi.cpp
target_include_directories(solve PRIVATE ${NETGEN_PYTHON_INCLUDE_DIRS})

if(NETGEN_USE_PYTHON)
add_library( ngslib SHARED ngspy.cpp )

set_target_properties(ngslib PROPERTIES SUFFIX ".so")
set_target_properties(ngslib PROPERTIES PREFIX "")
target_link_libraries(ngslib PRIVATE "$<BUILD_INTERFACE:netgen_python>")
target_link_libraries(ngslib ${LAPACK_CMAKE_LINK_INTERFACE} "$<BUILD_INTERFACE:ngs_lapack>")

target_link_libraries(ngslib PUBLIC ngsolve)

if(WIN32)
    set_target_properties( ngslib PROPERTIES SUFFIX ".pyd" )
endif(WIN32)

set_target_properties(ngslib PROPERTIES INSTALL_RPATH "${NETGEN_RPATH_TOKEN}/../${NETGEN_PYTHON_RPATH}")
install(TARGETS ngslib DESTINATION ${NGSOLVE_INSTALL_DIR_PYTHON}/ngsolve COMPONENT ngsolve)
endif(NETGEN_USE_PYTHON)

target_link_libraries( solve PUBLIC ngcomp ngstd ngfem ngla ngbla nglib ${MPI_CXX_LIBRARIES} PRIVATE "$<BUILD_INTERFACE:netgen_python>" ${VTUNE_LIB} ${PARDISO_LIBS})
target_link_libraries( solve ${LAPACK_CMAKE_LINK_INTERFACE} "$<BUILD_INTERFACE:ngs_lapack>")
install( TARGETS solve ${ngs_install_dir} )

add_library( ngsolve SHARED ${ngsolve_object_libs} ngsolve.cpp)
target_link_libraries(ngsolve PUBLIC nglib ngcore solve ngcomp ngfem ngbla ngla ngstd PRIVATE "$<BUILD_INTERFACE:netgen_python>" )
install( TARGETS ngsolve ${ngs_install_dir} )

if (WIN32)
    target_link_libraries(ngsolve PUBLIC ${MPI_CXX_LIBRARIES} PRIVATE "$<BUILD_INTERFACE:netgen_python>" ${PARDISO_LIBS} ${UMFPACK_LIBRARIES} )
    target_link_libraries(ngsolve ${LAPACK_CMAKE_LINK_INTERFACE} "$<BUILD_INTERFACE:ngs_lapack>")
endif()

if(WIN32)
    set_target_properties( ngsolve PROPERTIES OUTPUT_NAME libngsolve )
    enable_language(RC)
    target_sources( ngsolve PRIVATE ../cmake/resources/ngsolve.rc )
    set_target_properties( ngsolve PROPERTIES RULE_LAUNCH_COMPILE "" )
endif(WIN32)

install( FILES solve.hpp DESTINATION ${NGSOLVE_INSTALL_DIR_INCLUDE} COMPONENT ngsolve_devel )
