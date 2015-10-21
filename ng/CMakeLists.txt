set(netgen_sources demoview.cpp ngappinit.cpp onetcl.cpp parallelfunc.cpp ngpkg.cpp ../libsrc/stlgeom/stlpkg.cpp ../libsrc/visualization/visualpkg.cpp ../libsrc/csg/csgpkg.cpp ../libsrc/geom2d/geom2dpkg.cpp)

if(USE_GUI)
    if(WIN32)
      # add icon to netgen executable
      enable_language(RC)
      set(netgen_sources ${netgen_sources} ../windows/netgen.rc)
    endif(WIN32)

    add_executable(netgen ${netgen_sources})
    target_link_libraries( netgen nglib ${ZLIB_LIBRARIES} ${LIBTOGL} ${TK_LIBRARY} ${TCL_LIBRARY} ${JPEG_LIBRARIES} ${FFMPEG_LIBRARIES} ${X11_X11_LIB})
    if(NOT WIN32)
      target_link_libraries( netgen mesh stlvis stl geom2dvis interface geom2d csg stl visual csgvis )
    endif(NOT WIN32)

    install(TARGETS netgen    ${ng_install_dir})

endif(USE_GUI)

if(USE_PYTHON)
    add_library(ngpy SHARED netgenpy.cpp)
    target_link_libraries( ngpy nglib )
    if(APPLE)
        set_target_properties( ngpy PROPERTIES SUFFIX ".so")
    elseif(WIN32)
        set_target_properties( ngpy PROPERTIES SUFFIX ".pyd")
        set_target_properties( ngpy PROPERTIES OUTPUT_NAME "libngpy")
    endif()
    install(TARGETS ngpy ${ng_install_dir})
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

