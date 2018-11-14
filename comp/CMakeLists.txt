include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})

if (NETGEN_USE_GUI)
    add_definitions(-DTCL)
endif (NETGEN_USE_GUI)

if(NOT WIN32)
    set(ngcomp_object_libs $<TARGET_OBJECTS:ngmg> )
endif(NOT WIN32)

add_library( ngcomp ${NGS_LIB_TYPE}
        ${ngcomp_object_libs}
        bilinearform.cpp facetfespace.cpp fespace.cpp 
        gridfunction.cpp h1hofespace.cpp hcurlhdivfes.cpp hcurlhofespace.cpp 
        hdivfes.cpp hdivhofespace.cpp hdivhosurfacefespace.cpp hierarchicalee.cpp l2hofespace.cpp     
        linearform.cpp meshaccess.cpp ngsobject.cpp postproc.cpp	     
        preconditioner.cpp vectorfacetfespace.cpp numberfespace.cpp bddc.cpp h1amg.cpp
        hypre_precond.cpp hdivdivfespace.cpp hdivdivsurfacespace.cpp hcurlcurlfespace.cpp tpfes.cpp 
        python_comp.cpp python_comp_mesh.cpp ../fem/python_fem.cpp basenumproc.cpp pde.cpp pdeparser.cpp vtkoutput.cpp
        periodic.cpp hypre_ams_precond.cpp facetsurffespace.cpp compressedfespace.cpp
        )

target_compile_definitions(ngcomp PUBLIC ${NGSOLVE_COMPILE_DEFINITIONS})
target_compile_definitions(ngcomp PRIVATE ${NGSOLVE_COMPILE_DEFINITIONS_PRIVATE})
target_compile_options(ngcomp PUBLIC ${NGSOLVE_COMPILE_OPTIONS})
target_include_directories(ngcomp PUBLIC ${NGSOLVE_INCLUDE_DIRS})
target_include_directories(ngcomp PRIVATE ${NETGEN_TCL_INCLUDE_PATH} PRIVATE ${NETGEN_PYTHON_INCLUDE_DIRS})

if(NOT WIN32)
    target_link_libraries (ngcomp PUBLIC interface ngfem ngla ngbla ngstd ${MPI_CXX_LIBRARIES} PRIVATE ${NETGEN_PYTHON_LIBRARIES} ${HYPRE_LIBRARIES})
    target_link_libraries(ngcomp ${LAPACK_CMAKE_LINK_INTERFACE} ${LAPACK_LIBRARIES})
    install( TARGETS ngcomp ${ngs_install_dir} )
endif(NOT WIN32)

install( FILES
        bilinearform.hpp comp.hpp facetfespace.hpp	   
        fespace.hpp gridfunction.hpp h1hofespace.hpp hcurlhdivfes.hpp	   
        hcurlhofespace.hpp hdivfes.hpp hdivhofespace.hpp hdivhosurfacefespace.hpp		   	   
        l2hofespace.hpp hdivdivsurfacespace.hpp tpfes.hpp linearform.hpp meshaccess.hpp ngsobject.hpp	   
        postproc.hpp preconditioner.hpp vectorfacetfespace.hpp hypre_precond.hpp 
        pde.hpp numproc.hpp vtkoutput.hpp pmltrafo.hpp periodic.hpp  hypre_ams_precond.hpp facetsurffespace.hpp compressedfespace.hpp
        DESTINATION ${NGSOLVE_INSTALL_DIR_INCLUDE}
        COMPONENT ngsolve_devel
       )


