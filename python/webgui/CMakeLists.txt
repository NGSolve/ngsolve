include (ExternalProject)

set(sdir ${CMAKE_CURRENT_SOURCE_DIR})
set(bdir ${CMAKE_CURRENT_SOURCE_DIR}/js/package)

if(BUILD_JUPYTER_WIDGETS)
  find_program(NPM NAMES npm)
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
      ${bdir}/standalone.js
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
      COMMAND ${NPM} run build
      COMMAND cmake -E remove -f "*.tgz"
      COMMAND ${NPM} pack
      COMMAND cmake -E tar xvf "*.tgz"
      COMMAND cmake -E touch_nocreate "${bdir}/standalone.js"
      COMMAND cmake -E touch_nocreate "${bdir}/standalone.js"
      COMMAND cmake -E touch_nocreate "${bdir}/nbextension/static/index.js.map"
      COMMAND cmake -E touch_nocreate "${bdir}/nbextension/static/index.js"
      COMMAND cmake -E touch_nocreate "${bdir}/labextension/plugin.js.map"
      COMMAND cmake -E touch_nocreate "${bdir}/labextension/plugin.js"

      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/js
      )

  else(BUILD_JUPYTER_WIDGETS)
    # download prebuild widgets library

    # read lib version from package.json
    execute_process(COMMAND ${NETGEN_PYTHON_EXECUTABLE} -c
      "import json;print(json.load(open('package.json'))['version'])"
      OUTPUT_VARIABLE WIDGETS_VERSION
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/js
      OUTPUT_STRIP_TRAILING_WHITESPACE
    )

    ExternalProject_Add(
      download_widgets
      URL "https://registry.npmjs.org/ngsolve_jupyter_widgets/-/ngsolve_jupyter_widgets-${WIDGETS_VERSION}.tgz"
      URL_HASH SHA1=d660331f82d7d81e5b0a820a92db532a5db13d95
      DOWNLOAD_NO_PROGRESS ON
      SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/js/package
      DOWNLOAD_DIR ${CMAKE_SOURCE_DIR}/external_dependencies/
      BUILD_IN_SOURCE 1
      CONFIGURE_COMMAND cmake -E make_directory "${bdir}"
      BUILD_COMMAND cmake -E touch_nocreate "${bdir}/standalone.js"
      INSTALL_COMMAND cmake -E copy_directory . "${bdir}"
      UPDATE_COMMAND ""
      TEST_COMMAND ""
      )

    # adds dependency to behave like the content was compiled
    add_custom_command(OUTPUT ${bdir}/standalone.js COMMAND ; DEPENDS download_widgets)
endif(BUILD_JUPYTER_WIDGETS)

add_custom_command(OUTPUT webgui.py
  COMMAND ${NETGEN_PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/build.py ${sdir} ${bdir}
  DEPENDS
    build.py
    webgui_template.py
    ${bdir}/standalone.js
 )

add_custom_target(js_widgets ALL DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/webgui.py)

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/webgui.py DESTINATION ${NGSOLVE_INSTALL_DIR_PYTHON}/ngsolve COMPONENT ngsolve)

# jupyter notebook extension
install(DIRECTORY ${bdir}/nbextension DESTINATION ${NGSOLVE_INSTALL_DIR_PYTHON}/ngsolve COMPONENT ngsolve)
install(FILES js/extension.js DESTINATION ${NGSOLVE_INSTALL_DIR_PYTHON}/ngsolve/nbextension/static COMPONENT ngsolve)

# jupyter lab extension
install(DIRECTORY ${bdir}/labextension DESTINATION ${NGSOLVE_INSTALL_DIR_PYTHON}/ngsolve/labextension COMPONENT ngsolve)
install(FILES js/package/package.json DESTINATION ${NGSOLVE_INSTALL_DIR_PYTHON}/ngsolve/labextension COMPONENT ngsolve)

find_program(JUPYTER_NBEXTENSION NAMES jupyter-nbextension)
if(JUPYTER_NBEXTENSION)
  install(CODE "execute_process(COMMAND ${JUPYTER_NBEXTENSION} install --user --py ngsolve)")
  install(CODE "execute_process(COMMAND ${JUPYTER_NBEXTENSION} enable --user --py ngsolve)")
  message(STATUS "Found jupyter nbextension: ${JUPYTER_NBEXTENSION} -- ngsolve extension will be installed")
endif(JUPYTER_NBEXTENSION)

