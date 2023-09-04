set(NGSOLVE_USE_MKL $<BOOL:${USE_MKL}>)
configure_file(config.py ${CMAKE_CURRENT_BINARY_DIR}/config.py @ONLY)

if(NETGEN_USE_PYTHON)
    install (FILES
            ${CMAKE_CURRENT_BINARY_DIR}/config.py
            __expr.py internal.py __console.py webgui.py
            __init__.py utils.py solvers.py eigenvalues.py meshes.py
            krylovspace.py nonlinearsolvers.py bvp.py timing.py TensorProductTools.py
            ngs2petsc.py ngscxx.py config_cli.py
            DESTINATION ${NGSOLVE_INSTALL_DIR_PYTHON}/ngsolve
            COMPONENT ngsolve
            )

# build stub files for pybind11 packages
if(BUILD_STUB_FILES)
execute_process(COMMAND ${NETGEN_PYTHON_EXECUTABLE} -c "from importlib.metadata import version; print(version('pybind11-stubgen'))" OUTPUT_VARIABLE stubgen_version RESULT_VARIABLE pybind11_stubgen)
if(pybind11_stubgen AND NOT ${pybind11_stubgen} EQUAL 0)
  message(WARNING "pybind11-stubgen not found, if you want to create stub files
for better autocompletion support install it with pip.")
else()
  if(stubgen_version LESS "1.0")
    message(WARNING "pybind11-stubgen version is too old, if you want to create stub files for better autocompletion support upgrade it with pip.")
  else()
  message("-- Found pybind11-stubgen version: ${stubgen_version}")
  install(CODE "execute_process(COMMAND ${NETGEN_PYTHON_EXECUTABLE} -m pybind11_stubgen --ignore-all-errors ngsolve)")
  install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../stubs/ngsolve/ DESTINATION ${NGSOLVE_INSTALL_DIR_PYTHON}/ngsolve/ COMPONENT ngsolve)
  endif()
endif()
endif(BUILD_STUB_FILES)

endif(NETGEN_USE_PYTHON)
