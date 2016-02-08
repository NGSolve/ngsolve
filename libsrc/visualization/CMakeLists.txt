add_definitions(-DNGINTERFACE_EXPORTS)
install(FILES soldata.hpp DESTINATION include COMPONENT netgen_devel )

# AM_CPPFLAGS =  $(MPI_INCLUDES) -I$(top_srcdir)/libsrc/include  -DOPENGL -D$(TOGL_WINDOWINGSYSTEM) $(OCCFLAGS) $(TCL_INCLUDES)

if(USE_GUI)
    set( LIB_VISUAL_SOURCES meshdoc.cpp mvdraw.cpp vsfieldlines.cpp vsmesh.cpp vssolution.cpp importsolution.cpp )
else(USE_GUI)
    set( LIB_VISUAL_SOURCES visual_dummy.cpp )
endif(USE_GUI)

add_library(visual ${NG_LIB_TYPE} ${LIB_VISUAL_SOURCES})

if(NOT WIN32)
    target_link_libraries( visual ${PYTHON_LIBS} ${MPI_CXX_LIBRARIES} ${OPENGL_LIBRARIES} ${TCL_LIBRARY} )
    install( TARGETS visual ${ng_install_dir})
endif(NOT WIN32)





