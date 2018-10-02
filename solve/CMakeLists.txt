include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})

if(NOT WIN32)
endif(NOT WIN32)

add_library( solve ${NGS_LIB_TYPE}
        bvp.cpp generalizedEVP.cpp numproc.cpp evp.cpp 
        hyperbolic.cpp numprocee.cpp python_solve.cpp python_mpi.cpp)

target_compile_definitions(solve PUBLIC ${NGSOLVE_COMPILE_DEFINITIONS})
target_compile_definitions(solve PRIVATE ${NGSOLVE_COMPILE_DEFINITIONS_PRIVATE})
target_compile_options(solve PUBLIC ${NGSOLVE_COMPILE_OPTIONS})
target_include_directories(solve PUBLIC ${NGSOLVE_INCLUDE_DIRS})

if(NETGEN_USE_PYTHON)
add_library( ngslib SHARED ngspy.cpp )

set_target_properties(ngslib PROPERTIES SUFFIX ".so")
set_target_properties(ngslib PROPERTIES PREFIX "")
target_link_libraries(ngslib PUBLIC ${NETGEN_PYTHON_LIBRARIES})
target_link_libraries(ngslib PRIVATE ngmg ${LAPACK_CMAKE_LINK_INTERFACE} ${LAPACK_LIBRARIES})

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

if(NOT WIN32)
    target_link_libraries( solve PUBLIC ngcomp ngstd ngfem ngla ngbla nglib visual ${MPI_CXX_LIBRARIES} ${NETGEN_PYTHON_LIBRARIES} ${VTUNE_LIB} ${PARDISO_LIB})
    target_link_libraries( solve ${LAPACK_CMAKE_LINK_INTERFACE} ${LAPACK_LIBRARIES})
    target_link_libraries( solve PRIVATE ngmg )
    install( TARGETS solve ${ngs_install_dir} )
endif(NOT WIN32)
# target_link_libraries(ngslib ngsolve)



if (NETGEN_USE_GUI)
    if(WIN32)
        set(ngsolve_object_libs
            $<TARGET_OBJECTS:solve> 
            $<TARGET_OBJECTS:ngcomp> 
            $<TARGET_OBJECTS:ngstd> 
            $<TARGET_OBJECTS:ngfem> 
            $<TARGET_OBJECTS:ngla> 
            $<TARGET_OBJECTS:ngbla> 
#             $<TARGET_OBJECTS:interface> 
            )
    endif(WIN32)
    add_library( ngsolve SHARED ${ngsolve_object_libs} ngsolve.cpp shapetester.cpp )

    target_link_libraries(ngsolve PUBLIC nglib ${MPI_CXX_LIBRARIES} ${NETGEN_PYTHON_LIBRARIES} PRIVATE ${PARDISO_LIB} ${UMFPACK_LIBRARIES} ngmg)
    target_link_libraries(ngsolve ${LAPACK_CMAKE_LINK_INTERFACE} ${LAPACK_LIBRARIES})
    target_compile_definitions(ngsolve PUBLIC ${NGSOLVE_COMPILE_DEFINITIONS})
    target_compile_definitions(ngsolve PRIVATE ${NGSOLVE_COMPILE_DEFINITIONS_PRIVATE})
    target_compile_options(ngsolve PUBLIC ${NGSOLVE_COMPILE_OPTIONS})
    target_include_directories(ngsolve PUBLIC ${NGSOLVE_INCLUDE_DIRS})
    target_include_directories(ngsolve PRIVATE ${NETGEN_TCL_INCLUDE_PATH})

    if(NOT WIN32)
        target_link_libraries(ngsolve PUBLIC solve ngcomp ngstd ngfem ngla ngbla)
    endif(NOT WIN32)
    install( TARGETS ngsolve ${ngs_install_dir} )
    if(NOT APPLE)
        target_link_libraries( ngsolve PUBLIC gui)
    endif()
    if(WIN32)
        set_target_properties( ngsolve PROPERTIES OUTPUT_NAME libngsolve )
    endif(WIN32)
endif (NETGEN_USE_GUI)


if(NOT WIN32)
    add_executable( ngs ngs.cpp)
    target_link_libraries(ngs solve)
    install( TARGETS ngs DESTINATION ${NGSOLVE_INSTALL_DIR_BIN} COMPONENT ngsolve )
endif(NOT WIN32)

install( FILES solve.hpp DESTINATION ${NGSOLVE_INSTALL_DIR_INCLUDE} COMPONENT ngsolve_devel )
