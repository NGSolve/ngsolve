include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})

if(NOT WIN32)
    set(solve_object_libs
        $<TARGET_OBJECTS:ngmg> 
        $<TARGET_OBJECTS:parallel> 
        )
    set(ngslib_object_libs
        $<TARGET_OBJECTS:ngmg> 
        $<TARGET_OBJECTS:parallel> 
        )
endif(NOT WIN32)

add_library( solve ${NGS_LIB_TYPE}
        bvp.cpp generalizedEVP.cpp numproc.cpp evp.cpp 
        hyperbolic.cpp numprocee.cpp python_solve.cpp python_mpi.cpp
        ${solve_object_libs}
        )

if(USE_PYTHON)
add_library( ngslib SHARED ngspy.cpp )

set_target_properties(ngslib PROPERTIES SUFFIX ".so")
set_target_properties(ngslib PROPERTIES PREFIX "")
target_link_libraries(ngslib ${PYTHON_LIBS} ${LAPACK_LIBRARIES})

if(USE_GUI)
    target_link_libraries(ngslib ngsolve)
else(USE_GUI)
    target_link_libraries(ngslib solve)
endif(USE_GUI)

if(WIN32)
    set_target_properties( ngslib PROPERTIES SUFFIX ".pyd" )
endif(WIN32)

install( TARGETS ngslib ${ngs_install_dir} )
endif(USE_PYTHON)

if(NOT WIN32)
    target_link_libraries( solve ngcomp ngstd ngfem ngla ngbla interface visual ${MPI_CXX_LIBRARIES} ${PYTHON_LIBS} ${LAPACK_LIBRARIES} ${VTUNE_LIB} ${PARDISO_LIB})
    install( TARGETS solve ${ngs_install_dir} )
endif(NOT WIN32)
# target_link_libraries(ngslib ngsolve)



if (USE_GUI)
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
    add_library( ngsolve SHARED
            $<TARGET_OBJECTS:ngmg>
            $<TARGET_OBJECTS:parallel> 
            ${ngsolve_object_libs}
            ngsolve.cpp shapetester.cpp 
            )
    target_link_libraries(ngsolve interface ${MPI_CXX_LIBRARIES} ${PYTHON_LIBS} ${LAPACK_LIBRARIES} ${TCL_LIBRARY} ${PARDISO_LIB} ${MKL_LIBRARIES} )
    if(NOT WIN32)
        target_link_libraries(ngsolve solve ngcomp ngstd ngfem ngla ngbla)
    endif(NOT WIN32)
    install( TARGETS ngsolve ${ngs_install_dir} )
    if(WIN32)
        set_target_properties( ngsolve PROPERTIES OUTPUT_NAME libngsolve )
    endif(WIN32)
endif (USE_GUI)


if(NOT WIN32 AND NOT APPLE)
    add_executable( ngs ngs.cpp)
    target_link_libraries(ngs solve interface ${TCL_LIBRARY} ${PYTHON_LIBS} ${LAPACK_LIBRARIES})
    install( TARGETS ngs DESTINATION bin COMPONENT ngsolve )
endif(NOT WIN32 AND NOT APPLE)

install( FILES solve.hpp DESTINATION include  COMPONENT ngsolve_devel )
