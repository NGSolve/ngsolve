include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})

if(NGS_GOLD_DIR)
    add_definitions(-DGOLD -DJS)
    include_directories(BEFORE ${NGS_GOLD_DIR})
endif(NGS_GOLD_DIR)



add_library(ngfem ${NGS_LIB_TYPE}
        bdbequations.cpp maxwellintegrator.cpp
        hdiv_equations.cpp h1hofe.cpp h1lofe.cpp l2hofe.cpp hcurlhofe.cpp
        hcurlfe.cpp vectorfacetfe.cpp hdivhofe.cpp recursive_pol_trig.cpp
        coefficient.cpp integrator.cpp specialelement.cpp elementtopology.cpp
        intrule.cpp fastmat.cpp finiteelement.cpp elementtransformation.cpp
        scalarfe.cpp generic_recpol.cpp hdivfe.cpp recursive_pol.cpp
        hybridDG.cpp diffop.cpp l2hofefo.cpp h1hofefo.cpp
        facethofe.cpp python_fem.cpp DGIntegrators.cpp pml.cpp
        h1hofe_segm.cpp h1hofe_trig.cpp  symbolicintegrator.cpp
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
    target_link_libraries(ngfem ${PYTHON_LIBS} ${LAPACK_LIBRARIES})
    install( TARGETS ngfem ${ngs_install_dir} )

    if(NOT USE_NATIVE_ARCH)
        add_library(ngfem_avx2 ${NGS_LIB_TYPE} $<TARGET_PROPERTY:ngfem,SOURCES> )
        target_link_libraries(ngfem_avx2 $<TARGET_PROPERTY:ngfem,INTERFACE_LINK_LIBRARIES> )
        target_compile_options(ngfem_avx2 PUBLIC ${AVX2_FLAG})
        install( TARGETS ngfem_avx2 ${ngs_install_dir} )
    endif(NOT USE_NATIVE_ARCH)

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
        diffop_impl.hpp hcurlhofe_impl.hpp thcurlfe.hpp
        thcurlfe_impl.hpp symbolicintegrator.hpp code_generation.hpp
        DESTINATION include
        COMPONENT ngsolve_devel
       )

# install (FILES ngfemcf.py DESTINATION python)
