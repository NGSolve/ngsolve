add_definitions(-DNGINTERFACE_EXPORTS)
add_library(interface ${NG_LIB_TYPE}
#         $<TARGET_OBJECTS:mesh>
#         $<TARGET_OBJECTS:csg>
#         $<TARGET_OBJECTS:visual>
#         $<TARGET_OBJECTS:la>
#         $<TARGET_OBJECTS:gprim>
#         $<TARGET_OBJECTS:gen>
        nginterface.cpp nginterface_v2.cpp 
        read_fnf_mesh.cpp readtetmesh.cpp readuser.cpp writeabaqus.cpp writediffpack.cpp 
        writedolfin.cpp writeelmer.cpp writefeap.cpp writefluent.cpp writegmsh.cpp writejcm.cpp 
        writepermas.cpp writetecplot.cpp writetet.cpp writetochnog.cpp writeuser.cpp
        wuchemnitz.cpp writegmsh2.cpp writeOpenFOAM15x.cpp 
        )

# TARGET_LINK_LIBRARIES(interface ${MPI_CXX_LIBRARIES} ${PYTHON_LIBS})
# TARGET_LINK_LIBRARIES( interface ${LIBPTHREAD} ${OCC_LIBRARIES} ${LIBTOGL} ${MPI_CXX_LIBRARIES} ${OPENGL_LIBRARIES} ${TK_LIBRARY} ${TCL_LIBRARY} ${CMAKE_THREAD_LIBS_INIT} ${X11_Xmu_LIB} ${JPEG_LIBRARIES} ${MKL_LIBRARIES} ${ZLIB_LIBRARIES} )
if(NOT WIN32)
    target_link_libraries(interface mesh csg geom2d)
    if(USE_GUI)
        target_link_libraries(interface visual)
    endif(USE_GUI)
    install( TARGETS interface ${ng_install_dir})
endif(NOT WIN32)

