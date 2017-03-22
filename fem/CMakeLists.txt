include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})

if(NGS_GOLD_DIR)
    add_definitions(-DGOLD -DJS)
    include_directories(BEFORE ${NGS_GOLD_DIR})
endif(NGS_GOLD_DIR)



add_library(ngfem ${NGS_LIB_TYPE}
        bdbequations.cpp maxwellintegrator.cpp
        hdiv_equations.cpp h1hofe.cpp h1lofe.cpp l2hofe.cpp  l2hofe_trig.cpp l2hofe_segm.cpp l2hofe_tet.cpp hcurlhofe.cpp
        hcurlhofe_hex.cpp hcurlhofe_tet.cpp hcurlhofe_prism.cpp hcurlhofe_pyramid.cpp
        hcurlfe.cpp vectorfacetfe.cpp hdivhofe.cpp recursive_pol_trig.cpp
        coefficient.cpp integrator.cpp specialelement.cpp elementtopology.cpp
        intrule.cpp fastmat.cpp finiteelement.cpp elementtransformation.cpp
        scalarfe.cpp generic_recpol.cpp hdivfe.cpp recursive_pol.cpp
        hybridDG.cpp diffop.cpp l2hofefo.cpp h1hofefo.cpp
        facethofe.cpp python_fem.cpp DGIntegrators.cpp pml.cpp
        h1hofe_segm.cpp h1hofe_trig.cpp symbolicintegrator.cpp tpdiffop.cpp
        tensorproductintegrator.cpp
        )

if(USE_CUDA)
    cuda_add_library( cuda_fem STATIC
            test.cu test1.cu fem_kernels.cu
            OPTIONS ${CUDA_OPTIONS}
            )
    target_link_libraries(ngfem cuda_fem)
    target_link_libraries(ngfem ${CUDA_LIBRARIES} ${CUDA_CUBLAS_LIBRARIES} ${CUDA_cusparse_LIBRARY})
endif(USE_CUDA)

if(NOT WIN32)
    target_link_libraries(ngfem ngbla ngstd ${NETGEN_PYTHON_LIBRARIES} ${LAPACK_LIBRARIES})
    install( TARGETS ngfem ${ngs_install_dir} )
endif(NOT WIN32)

install( FILES
        bdbequations.hpp hcurl_equations.hpp		
        elasticity_equations.hpp diffop.hpp bdbintegrator.hpp coefficient.hpp
        elementtopology.hpp elementtransformation.hpp facetfe.hpp	
        facethofe.hpp fastmat.hpp fem.hpp finiteelement.hpp generic_recpol.hpp	
        h1hofefo.hpp h1hofefo_impl.hpp h1hofe.hpp h1lofe.hpp hcurlfe.hpp
        hcurlhofe.hpp hcurllofe.hpp hdiv_equations.hpp hdivfe.hpp hdivhofe.hpp
        integrator.hpp intrule.hpp l2hofefo.hpp l2hofe.hpp recursive_pol.hpp
        recursive_pol_tet.hpp recursive_pol_trig.hpp scalarfe.hpp	
        specialelement.hpp thdivfe.hpp tscalarfe.hpp vectorfacetfe.hpp	
        hdivlofe.hpp hdivhofefo.hpp pml.hpp precomp.hpp h1hofe_impl.hpp	
        hdivhofe_impl.hpp tscalarfe_impl.hpp thdivfe_impl.hpp l2hofe_impl.hpp
        diffop_impl.hpp hcurlhofe_impl.hpp thcurlfe.hpp tpdiffop.hpp tpintrule.hpp
        thcurlfe_impl.hpp symbolicintegrator.hpp code_generation.hpp 
        tensorproductintegrator.hpp fe_interfaces.hpp
        DESTINATION include
        COMPONENT ngsolve_devel
       )

# install (FILES ngfemcf.py DESTINATION python)
