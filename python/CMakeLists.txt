if(NETGEN_USE_PYTHON)
    install (FILES
            __expr.py internal.py __console.py
            __init__.py utils.py solvers.py eigenvalues.py meshes.py
            krylovspace.py nonlinearsolvers.py bvp.py timing.py TensorProductTools.py
            DESTINATION ${NGSOLVE_INSTALL_DIR_PYTHON}/ngsolve
            COMPONENT ngsolve
            )

# build stub files for pybind11 packages
if(BUILD_STUB_FILES)
find_program(PYBIND11_STUBS NAMES pybind11-stubgen)
if(PYBIND11_STUBS)
  message("-- Found pybind11-stubgen: ${PYBIND11_STUBS}")
  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../stubs/ngsolve-stubs/)
  install(CODE "execute_process(COMMAND ${PYBIND11_STUBS} --no-setup-py ngsolve)")
  install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../stubs/ngsolve-stubs/ DESTINATION ${NGSOLVE_INSTALL_DIR_PYTHON}/ngsolve/ COMPONENT ngsolve)
else(PYBIND11_STUBS)
  message(WARNING "pybind11-stubgen not found, if you want to create stub files
for better autocompletion support install it with pip.")
endif(PYBIND11_STUBS)
endif(BUILD_STUB_FILES)

add_subdirectory(webgui)

endif(NETGEN_USE_PYTHON)
