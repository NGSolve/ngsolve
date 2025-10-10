include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})

add_library( ngcomp ${NGS_LIB_TYPE}
        ${ngcomp_object_libs}
        bilinearform.cpp facetfespace.cpp fespace.cpp 
        gridfunction.cpp h1hofespace.cpp hcurlhdivfes.cpp hcurlhofespace.cpp 
        hdivfes.cpp hdivhofespace.cpp hdivhosurfacefespace.cpp l2hofespace.cpp     
        linearform.cpp meshaccess.cpp ngsobject.cpp postproc.cpp	     
        preconditioner.cpp tangentialfacetfespace.cpp
        normalfacetfespace.cpp normalfacetsurfacefespace.cpp
        numberfespace.cpp irspace.cpp h1lumping.cpp bddc.cpp h1amg.cpp
        hypre_precond.cpp hdivdivfespace.cpp hdivdivsurfacespace.cpp hcurlcurlfespace.cpp tpfes.cpp hcurldivfespace.cpp fesconvert.cpp
        python_comp.cpp python_comp_mesh.cpp ../fem/python_fem.cpp vtkoutput.cpp
        periodic.cpp discontinuous.cpp hidden.cpp reorderedfespace.cpp
        hypre_ams_precond.cpp facetsurffespace.cpp
        compressedfespace.cpp plateaufespace.cpp 
        globalinterfacespace.cpp globalspace.cpp
        ../multigrid/mgpre.cpp ../multigrid/prolongation.cpp
        ../multigrid/smoother.cpp contact.cpp localsolve.cpp interpolate.cpp
        webgui.cpp hcurlamg.cpp matrixfree.cpp
        )

target_include_directories(ngcomp PRIVATE ${NETGEN_PYTHON_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}/../ngstd ${CMAKE_CURRENT_SOURCE_DIR}/../linalg)

target_link_libraries (ngcomp PUBLIC ngfem ngla ngbla ngstd ${MPI_CXX_LIBRARIES} PRIVATE "$<BUILD_INTERFACE:netgen_python>" ${HYPRE_LIBRARIES})
target_link_libraries(ngcomp ${LAPACK_CMAKE_LINK_INTERFACE} "$<BUILD_INTERFACE:ngs_lapack>")
target_link_libraries(ngcomp PUBLIC netgen_lib)
install( TARGETS ngcomp ${ngs_install_dir} )

# disable contact pair visualization due to linking issues with netgen gui lib
# if(NETGEN_USE_GUI)
#   target_link_libraries(ngcomp PUBLIC nggui)
# endif(NETGEN_USE_GUI)

install( FILES
        bilinearform.hpp comp.hpp facetfespace.hpp	   
        fespace.hpp gridfunction.hpp h1hofespace.hpp hcurlhdivfes.hpp	   
        hcurlhofespace.hpp hdivfes.hpp hdivhofespace.hpp hdivhosurfacefespace.hpp		   	   
        l2hofespace.hpp hdivdivsurfacespace.hpp tpfes.hpp linearform.hpp meshaccess.hpp ngsobject.hpp	   
        postproc.hpp preconditioner.hpp tangentialfacetfespace.hpp
        normalfacetfespace.hpp normalfacetsurfacefespace.hpp hypre_precond.hpp h1amg.hpp
        irspace.hpp h1lumping.hpp vtkoutput.hpp pmltrafo.hpp periodic.hpp
        discontinuous.hpp hidden.hpp reorderedfespace.hpp
        hypre_ams_precond.hpp facetsurffespace.hpp
        compressedfespace.hpp plateaufespace.hpp globalinterfacespace.hpp globalspace.hpp
        python_comp.hpp fesconvert.hpp contact.hpp interpolate.hpp
        webgui.hpp hcurlcurlfespace.hpp
        DESTINATION ${NGSOLVE_INSTALL_DIR_INCLUDE}
        COMPONENT ngsolve_devel
       )


