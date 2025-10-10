include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})

add_library( ngsolve ${NGSOLVE_LIB_TYPE} ${ngsolve_object_libs} ngsolve.cpp python_solve.cpp)
target_include_directories(ngsolve PRIVATE ${NETGEN_PYTHON_INCLUDE_DIRS})
target_link_libraries(ngsolve PUBLIC netgen_lib netgen_core_lib solve ngcomp ngfem ngsbem ngbla ngla ngstd PRIVATE "$<BUILD_INTERFACE:netgen_python>" )
install( TARGETS ngsolve ${ngs_install_dir} )
add_library(solve ALIAS ngsolve)

if (WIN32)
    target_link_libraries(ngsolve PUBLIC ${MPI_CXX_LIBRARIES} PRIVATE ${VTUNE_LIB} ${PARDISO_LIBS} ${UMFPACK_LIBRARIES} )
    target_link_libraries(ngsolve ${LAPACK_CMAKE_LINK_INTERFACE} "$<BUILD_INTERFACE:ngs_lapack>")

    set_target_properties( ngsolve PROPERTIES OUTPUT_NAME libngsolve )
    enable_language(RC)
    target_sources( ngsolve PRIVATE ../cmake/resources/ngsolve.rc )
    set_target_properties( ngsolve PROPERTIES RULE_LAUNCH_COMPILE "" )
endif(WIN32)

if(NETGEN_USE_PYTHON)
    add_library( ngslib SHARED ngspy.cpp )

    set_target_properties(ngslib PROPERTIES SUFFIX ".so")
    set_target_properties(ngslib PROPERTIES PREFIX "")
    target_link_libraries(ngslib PRIVATE "$<BUILD_INTERFACE:netgen_python>")
    target_link_libraries(ngslib ${LAPACK_CMAKE_LINK_INTERFACE} "$<BUILD_INTERFACE:ngs_lapack>")

    target_link_libraries(ngslib PUBLIC ngsolve)

    if(WIN32)
        set_target_properties( ngslib PROPERTIES SUFFIX ".pyd" )
    else(WIN32)
        target_link_libraries(ngslib PRIVATE ngstd ngbla ngla ngfem ngcomp ngsbem)
    endif(WIN32)

    set_target_properties(ngslib PROPERTIES INSTALL_RPATH "${NETGEN_RPATH_TOKEN}/../${NETGEN_PYTHON_RPATH}")
    install(TARGETS ngslib DESTINATION ${NGSOLVE_INSTALL_DIR_PYTHON}/ngsolve COMPONENT ngsolve)
endif(NETGEN_USE_PYTHON)

install( FILES solve.hpp DESTINATION ${NGSOLVE_INSTALL_DIR_INCLUDE} COMPONENT ngsolve_devel )
