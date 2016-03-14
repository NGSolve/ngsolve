add_definitions(-DNGLIB_EXPORTS)

if(WIN32)
    set(nglib_objects 
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
          $<TARGET_OBJECTS:occ>
       )
endif(WIN32)

add_library(nglib SHARED nglib.cpp ${nglib_objects})
if(NOT WIN32)
      target_link_libraries( nglib mesh stlvis stl geom2dvis interface geom2d csg stl visual csgvis )
endif(NOT WIN32)

target_link_libraries( nglib ${OCC_LIBRARIES} ${MPI_CXX_LIBRARIES} ${OPENGL_LIBRARIES} ${TK_LIBRARY} ${TCL_LIBRARY} ${CMAKE_THREAD_LIBS_INIT} ${X11_Xmu_LIB} ${JPEG_LIBRARIES} ${MKL_LIBRARIES} ${ZLIB_LIBRARIES} ${OCC_LIBRARIES} )

if(USE_OCC AND NOT WIN32)
    target_link_libraries(nglib occ)
endif(USE_OCC AND NOT WIN32)

if(USE_PYTHON)
    target_link_libraries(nglib ${PYTHON_LIBS})
endif(USE_PYTHON)

# if(NOT WIN32)
#     add_executable(ng_vol ng_vol.cpp ${nglib_objects})
#     target_link_libraries( ng_vol nglib ${OCC_LIBRARIES} ${MPI_CXX_LIBRARIES} )
#     target_link_libraries( ng_vol nglib interface geom2d csg stl occ mesh )
# 
#     add_executable(ng_stl ng_stl.cpp ${nglib_objects})
#     target_link_libraries( ng_stl nglib ${OCC_LIBRARIES} ${OCC_OCAF_LIBRARIES} ${MPI_CXX_LIBRARIES} )
#     target_link_libraries( ng_stl nglib interface geom2d csg stl occ mesh )
# 
#     install(TARGETS nglib ng_vol ng_stl ${ng_install_dir})
# endif(NOT WIN32)

install(TARGETS nglib ${ng_install_dir})
