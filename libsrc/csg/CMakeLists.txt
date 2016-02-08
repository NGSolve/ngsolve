add_definitions(-DNGINTERFACE_EXPORTS)
add_library(csg ${NG_LIB_TYPE}
        algprim.cpp brick.cpp   
        bspline2d.cpp csgeom.cpp csgparser.cpp curve2d.cpp edgeflw.cpp	       
        explicitcurve2d.cpp extrusion.cpp gencyl.cpp genmesh.cpp identify.cpp  
        manifold.cpp meshsurf.cpp polyhedra.cpp revolution.cpp singularref.cpp 
        solid.cpp specpoin.cpp spline3d.cpp surface.cpp triapprox.cpp zrefine.cpp 
        python_csg.cpp
        )
if(APPLE)
    set_target_properties( csg PROPERTIES SUFFIX ".so")
endif(APPLE)

if(NOT WIN32)
    target_link_libraries(csg mesh ${PYTHON_LIBS})
    target_link_libraries(csg ${PYTHON_LIBS})
    install( TARGETS csg ${ng_install_dir})
endif(NOT WIN32)


if(USE_GUI)
    add_library(csgvis ${NG_LIB_TYPE} vscsg.cpp )
    if(NOT WIN32)
        target_link_libraries(csgvis csg visual)
        if(APPLE)
            set_target_properties( csgvis PROPERTIES SUFFIX ".so")
        endif(APPLE)
        install( TARGETS csgvis ${ng_install_dir})
    endif(NOT WIN32)
endif(USE_GUI)

