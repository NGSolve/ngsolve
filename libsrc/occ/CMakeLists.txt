add_library(occ ${NG_LIB_TYPE}
        Partition_Inter2d.cxx Partition_Inter3d.cxx 
        Partition_Loop.cxx Partition_Loop2d.cxx Partition_Loop3d.cxx Partition_Spliter.cxx 
        occconstruction.cpp occgenmesh.cpp occgeom.cpp occmeshsurf.cpp
        )

add_library(occvis ${NG_LIB_TYPE} vsocc.cpp)

if(NOT WIN32)
    target_link_libraries( occ ${OCC_LIBRARIES} ${PYTHON_LIBS})
    target_link_libraries( occvis occ )
    install( TARGETS occ occvis ${ng_install_dir})
endif(NOT WIN32)

