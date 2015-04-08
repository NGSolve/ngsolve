add_library(stl ${NG_LIB_TYPE} 
        meshstlsurface.cpp stlgeom.cpp stlgeomchart.cpp 
        stlgeommesh.cpp stlline.cpp stltool.cpp stltopology.cpp
        )

if(NOT WIN32)
    target_link_libraries( stl mesh  ${PYTHON_LIBS})
    target_link_libraries( stl ${PYTHON_LIBS})
    install( TARGETS stl DESTINATION lib COMPONENT netgen )
endif(NOT WIN32)

if(USE_GUI)
    add_library(stlvis ${NG_LIB_TYPE}
            stlpkg.cpp vsstl.cpp
            )
    if(NOT WIN32)
        target_link_libraries( stlvis stl )
        install( TARGETS stlvis DESTINATION lib COMPONENT netgen )
    endif(NOT WIN32)
endif(USE_GUI)

