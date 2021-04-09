if(NETGEN_USE_PYTHON)
    install (FILES
            __expr.py internal.py __console.py webgui.py
            __init__.py utils.py solvers.py eigenvalues.py meshes.py
            krylovspace.py nonlinearsolvers.py bvp.py timing.py TensorProductTools.py
            ngs2petsc.py
            DESTINATION ${NGSOLVE_INSTALL_DIR_PYTHON}/ngsolve
            COMPONENT ngsolve
            )

# build stub files for pybind11 packages
if(BUILD_STUB_FILES)
execute_process(COMMAND ${NETGEN_PYTHON_EXECUTABLE} -c "import pybind11_stubgen; print(pybind11_stubgen.__file__)" OUTPUT_VARIABLE stubgen_path RESULT_VARIABLE pybind11_stubgen)
if(pybind11_stubgen AND NOT ${pybind11_stubgen} EQUAL 0)
  message(WARNING "pybind11-stubgen not found, if you want to create stub files
for better autocompletion support install it with pip.")
else()
  message("-- Found pybind11-stubgen: ${stubgen_path}")
  install(CODE "execute_process(COMMAND ${NETGEN_PYTHON_EXECUTABLE} -m pybind11_stubgen --no-setup-py ngsolve)") 
  install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../stubs/ngsolve-stubs/ DESTINATION ${NGSOLVE_INSTALL_DIR_PYTHON}/ngsolve/ COMPONENT ngsolve)
endif()
endif(BUILD_STUB_FILES)

endif(NETGEN_USE_PYTHON)
