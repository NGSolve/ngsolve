include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})

add_library( solve ${NGS_LIB_TYPE}
        bvp.cpp generalizedEVP.cpp numproc.cpp evp.cpp 
        hyperbolic.cpp numprocee.cpp python_solve.cpp
        )
#  python_mpi.cpp
target_include_directories(solve PRIVATE ${NETGEN_PYTHON_INCLUDE_DIRS})

if(NETGEN_USE_PYTHON)
add_library( ngslib SHARED ngspy.cpp )

set_target_properties(ngslib PROPERTIES SUFFIX ".so")
set_target_properties(ngslib PROPERTIES PREFIX "")
target_link_libraries(ngslib PRIVATE netgen_python)
target_link_libraries(ngslib ${LAPACK_CMAKE_LINK_INTERFACE} ${LAPACK_LIBRARIES})

if(NETGEN_USE_GUI AND WIN32)
    target_link_libraries(ngslib PUBLIC ngsolve)
else(NETGEN_USE_GUI AND WIN32)
    target_link_libraries(ngslib PUBLIC solve)
endif(NETGEN_USE_GUI AND WIN32)

if(WIN32)
    set_target_properties( ngslib PROPERTIES SUFFIX ".pyd" )
endif(WIN32)

set_target_properties(ngslib PROPERTIES INSTALL_RPATH "${NETGEN_RPATH_TOKEN}/../${NETGEN_PYTHON_RPATH}")
install(TARGETS ngslib DESTINATION ${NGSOLVE_INSTALL_DIR_PYTHON}/ngsolve COMPONENT ngsolve)
endif(NETGEN_USE_PYTHON)

target_link_libraries( solve PUBLIC ngcomp ngstd ngfem ngla ngbla nglib visual ${MPI_CXX_LIBRARIES} PRIVATE netgen_python ${VTUNE_LIB} ${PARDISO_LIBS})
target_link_libraries( solve ${LAPACK_CMAKE_LINK_INTERFACE} ${LAPACK_LIBRARIES})
install( TARGETS solve ${ngs_install_dir} )

if (NETGEN_USE_GUI)
    add_library( ngsolve SHARED
            ${ngsolve_object_libs}
            ngsolve.cpp shapetester.cpp 
            )

    target_link_libraries(ngsolve PUBLIC nglib gui ${MPI_CXX_LIBRARIES} PRIVATE netgen_python ${PARDISO_LIBS} ${UMFPACK_LIBRARIES} )
    target_link_libraries(ngsolve ${LAPACK_CMAKE_LINK_INTERFACE} ${LAPACK_LIBRARIES})
    target_include_directories(ngsolve PRIVATE ${NETGEN_TCL_INCLUDE_PATH})

    target_link_libraries(ngsolve PUBLIC solve ngcomp ngstd ngfem ngla ngbla)
    install( TARGETS ngsolve ${ngs_install_dir} )

    if(NOT APPLE)
        target_link_libraries( ngsolve PRIVATE ${NETGEN_TCL_LIBRARY} ${NETGEN_TK_LIBRARY} )
    endif(NOT APPLE)
    if(WIN32)
        set_target_properties( ngsolve PROPERTIES OUTPUT_NAME libngsolve )
    endif(WIN32)
endif (NETGEN_USE_GUI)


if(NOT WIN32)
    add_executable( ngs ngs.cpp)
    target_link_libraries(ngs PUBLIC solve netgen_python ${NETGEN_PYTHON_LIBRARIES})
    install( TARGETS ngs DESTINATION ${NGSOLVE_INSTALL_DIR_BIN} COMPONENT ngsolve )
endif(NOT WIN32)

install( FILES solve.hpp DESTINATION ${NGSOLVE_INSTALL_DIR_INCLUDE} COMPONENT ngsolve_devel )
