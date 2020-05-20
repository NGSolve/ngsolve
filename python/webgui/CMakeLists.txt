set(sdir ${CMAKE_CURRENT_SOURCE_DIR})
set(bdir ${CMAKE_CURRENT_SOURCE_DIR}/js/package)
message("dirs ${sdir} ${bdir}")

if(BUILD_JUPYTER_WIDGETS)
  add_custom_command(OUTPUT ${sdir}/js/src/shaders.ts
    COMMAND ${NETGEN_PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/build_shaders.py ${CMAKE_CURRENT_SOURCE_DIR}
    DEPENDS
      build_shaders.py
      shader/clipping_vectors.vert
      shader/function.frag
      shader/splines.frag
      shader/splines.vert
      shader/trigsplines.vert
      shader/utils.h
      shader/vector_function.vert
    )

  add_custom_command(
    OUTPUT 
      ${bdir}/dist/index.js
      ${bdir}/dist/index.js.map
      ${bdir}/docs/source/_static/embed-bundle.js
      ${bdir}/docs/source/_static/embed-bundle.js.map
      ${bdir}/nbextension/static
      ${bdir}/nbextension/static/index.js
      ${bdir}/nbextension/static/index.js.map
      ${bdir}/standalone.js
      ${bdir}/standalone.js.map
    COMMAND npm install
    COMMAND npm run build
    DEPENDS
      js/package.json
      js/tsconfig.json
      js/webpack.config.js
      js/src/extension.ts
      js/src/index.ts
      js/src/plugin.ts
      js/src/scene.ts
      js/src/shaders.ts
      js/src/version.ts
      js/src/widget.ts
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/js
    )

else(BUILD_JUPYTER_WIDGETS)
  # download prebuild widgets library
endif(BUILD_JUPYTER_WIDGETS)

add_custom_command(OUTPUT webgui.py
  COMMAND ${NETGEN_PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/build.py ${sdir} ${bdir}
  DEPENDS
    build.py
    webgui_template.py
    ${bdir}/standalone.js
 )

add_custom_target(js_widgets ALL DEPENDS webgui.py)

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/webgui.py DESTINATION ${NGSOLVE_INSTALL_DIR_PYTHON}/ngsolve COMPONENT ngsolve)

# jupyter notebook extension
install(DIRECTORY ${bdir}/nbextension DESTINATION ${NGSOLVE_INSTALL_DIR_PYTHON}/ngsolve COMPONENT ngsolve)
install(FILES js/extension.js DESTINATION ${NGSOLVE_INSTALL_DIR_PYTHON}/ngsolve/nbextension/static COMPONENT ngsolve)

# jupyter lab extension
install(DIRECTORY ${sdir}/js/build/lib DESTINATION ${NGSOLVE_INSTALL_DIR_PYTHON}/ngsolve/labextension COMPONENT ngsolve)
install(FILES js/package.json DESTINATION ${NGSOLVE_INSTALL_DIR_PYTHON}/ngsolve/labextension COMPONENT ngsolve)
