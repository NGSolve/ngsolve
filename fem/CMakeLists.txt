include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})

if(NGS_GOLD_DIR)
    add_definitions(-DGOLD -DJS)
    include_directories(BEFORE ${NGS_GOLD_DIR})
endif(NGS_GOLD_DIR)



add_library(ngfem ${NGS_LIB_TYPE}
        bdbequations.cpp diffop_grad.cpp diffop_hesse.cpp
        diffop_id.cpp maxwellintegrator.cpp
        hdiv_equations.cpp h1hofe.cpp nodalhofe.cpp h1lofe.cpp l2hofe.cpp
        l2hofe_trig.cpp l2hofe_segm.cpp l2hofe_tet.cpp l2hofetp.cpp hcurlhofe.cpp
        hcurlhofe_hex.cpp hcurlhofe_tet.cpp hcurlhofe_prism.cpp hcurlhofe_pyramid.cpp
        hcurlfe.cpp vectorfacetfe.cpp normalfacetfe.cpp hdivhofe.cpp recursive_pol_trig.cpp
        coefficient.cpp coefficient_geo.cpp coefficient_stdmath.cpp coefficient_impl.hpp
        integrator.cpp specialelement.cpp elementtopology.cpp
        intrule.cpp finiteelement.cpp elementtransformation.cpp
        scalarfe.cpp hdivfe.cpp recursive_pol.cpp
        hybridDG.cpp diffop.cpp l2hofefo.cpp h1hofefo.cpp
        facethofe.cpp DGIntegrators.cpp pml.cpp
        h1hofe_segm.cpp h1hofe_trig.cpp h1hofe_quad.cpp h1hofe_tet.cpp h1hofe_prism.cpp h1hofe_pyramid.cpp h1hofe_hex.cpp
        hdivdivfe.cpp hcurlcurlfe.cpp symbolicintegrator.cpp tpdiffop.cpp
        newtonCF.cpp tensorproductintegrator.cpp code_generation.cpp
        voxelcoefficientfunction.cpp
        tensorcoefficient.cpp tensorcoefficient.hpp)

target_link_libraries(ngfem PUBLIC ngbla ngstd PRIVATE "$<BUILD_INTERFACE:netgen_python>")
target_link_libraries(ngfem ${LAPACK_CMAKE_LINK_INTERFACE} "$<BUILD_INTERFACE:ngs_lapack>")
install( TARGETS ngfem ${ngs_install_dir} )

install( FILES
        bdbequations.hpp hcurl_equations.hpp		
        elasticity_equations.hpp diffop.hpp bdbintegrator.hpp
        coefficient.hpp coefficient_impl.hpp coefficient_stdmath.hpp
        elementtopology.hpp elementtransformation.hpp facetfe.hpp	
        facethofe.hpp fem.hpp finiteelement.hpp
        shapefunction_utils.hpp hcurlfe_utils.hpp hdivfe_utils.hpp
        h1hofefo.hpp h1hofefo_impl.hpp h1hofe.hpp nodalhofe.hpp nodalhofe_impl.hpp h1lofe.hpp hcurlfe.hpp
        hcurlhofe.hpp hcurllofe.hpp hdivdivfe.hpp hdiv_equations.hpp hdivfe.hpp hdivhofe.hpp
        integrator.hpp integratorcf.hpp intrule.hpp l2hofefo.hpp l2hofe.hpp recursive_pol.hpp
        recursive_pol_tet.hpp recursive_pol_trig.hpp scalarfe.hpp	
        specialelement.hpp thdivfe.hpp tscalarfe.hpp vectorfacetfe.hpp normalfacetfe.hpp
        hdivlofe.hpp hdivhofefo.hpp pml.hpp precomp.hpp h1hofe_impl.hpp	
        hdivhofe_impl.hpp tscalarfe_impl.hpp thdivfe_impl.hpp
        l2hofe_impl.hpp hcurlcurlfe.hpp
        diffop_impl.hpp hcurlhofe_impl.hpp thcurlfe.hpp tpdiffop.hpp tpintrule.hpp
        thcurlfe_impl.hpp symbolicintegrator.hpp hcurlhdiv_dshape.hpp code_generation.hpp 
        tensorcoefficient.hpp tensorproductintegrator.hpp fe_interfaces.hpp python_fem.hpp
        voxelcoefficientfunction.hpp
        DESTINATION ${NGSOLVE_INSTALL_DIR_INCLUDE}
        COMPONENT ngsolve_devel
       )
