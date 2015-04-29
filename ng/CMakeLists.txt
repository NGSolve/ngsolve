add_definitions(-DNGINTERFACE_EXPORTS)

set(netgen_sources main.cpp )
set(netgenlib_sources demoview.cpp ngappinit.cpp onetcl.cpp parallelfunc.cpp ngpkg.cpp )

if(USE_GUI)
    if(WIN32)
    set(netgenlib_objects
            $<TARGET_OBJECTS:mesh>
            $<TARGET_OBJECTS:stlvis>
            $<TARGET_OBJECTS:stl>
            $<TARGET_OBJECTS:geom2dvis>
            $<TARGET_OBJECTS:interface>
            $<TARGET_OBJECTS:geom2d>
            $<TARGET_OBJECTS:csg>
            $<TARGET_OBJECTS:stl>
            $<TARGET_OBJECTS:gen>
            $<TARGET_OBJECTS:la>
            $<TARGET_OBJECTS:gprim>

            $<TARGET_OBJECTS:visual>
            $<TARGET_OBJECTS:csgvis>

#         $<TARGET_OBJECTS:occ>
#         $<TARGET_OBJECTS:occvis>
       )
    endif(WIN32)

    add_executable(netgen ${netgen_sources})
    add_library(netgenlib SHARED ${netgenlib_sources} ${netgenlib_objects})

    if(WIN32)
        set_target_properties( netgenlib PROPERTIES OUTPUT_NAME interface ) 
    else(WIN32)
        target_link_libraries( netgenlib visual csgvis csg interface mesh occ occvis togl)
    endif(WIN32)

    target_link_libraries( netgenlib ${LIBPTHREAD} ${OCC_LIBRARIES} ${LIBTOGL} ${MPI_CXX_LIBRARIES} ${OPENGL_LIBRARIES} ${TK_LIBRARY} ${TCL_LIBRARY} ${CMAKE_THREAD_LIBS_INIT} ${X11_Xmu_LIB} ${JPEG_LIBRARIES} ${MKL_LIBRARIES} ${ZLIB_LIBRARIES} ${PYTHON_LIBS} ${FFMPEG_LIBRARIES})
    target_link_libraries( netgen netgenlib )

    install(TARGETS netgenlib ${ng_install_dir})
    install(TARGETS netgen    ${ng_install_dir})

endif(USE_GUI)

if(USE_PYTHON)
    add_library(pynglib MODULE netgenpy.cpp)
    target_link_libraries( pynglib nglib )
    if(APPLE)
        set_target_properties( pynglib PROPERTIES SUFFIX ".so")
    elseif(WIN32)
        set_target_properties( pynglib PROPERTIES SUFFIX ".pyd")
    endif()
    set_target_properties( pynglib PROPERTIES OUTPUT_NAME "nglib")
    install(TARGETS pynglib ${ng_install_dir})
endif(USE_PYTHON)

install(FILES
    dialog.tcl menustat.tcl ngicon.tcl ng.tcl	
    ngvisual.tcl sockets.tcl drawing.tcl nghelp.tcl ngshell.tcl	  
    ngtesting.tcl parameters.tcl variables.tcl csgeom.tcl stlgeom.tcl 
    occgeom.tcl acisgeom.tcl netgen.ocf
    DESTINATION bin COMPONENT netgen)

if(USE_GUI AND NOT WIN32)
    add_subdirectory(Togl-1.7)
endif(USE_GUI AND NOT WIN32)

