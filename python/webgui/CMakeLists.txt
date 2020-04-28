add_custom_target(webgui.py ALL
  COMMAND ${NETGEN_PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/build.py ${CMAKE_CURRENT_SOURCE_DIR}
  SOURCES
    webgui.py
    render.js
    shader/clipping_vectors.vert
    shader/function.frag
    shader/splines.frag
    shader/splines.vert
    shader/trigsplines.vert
    shader/utils.h
    shader/vector_function.vert
  )

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/webgui.py DESTINATION ${NGSOLVE_INSTALL_DIR_PYTHON}/ngsolve COMPONENT ngsolve)
