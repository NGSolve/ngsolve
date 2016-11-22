install(DIRECTORY netgen/external_dependencies/pybind11/include/pybind11 DESTINATION include COMPONENT ngsolve_devel)
install(FILES netgen/external_dependencies/pybind11/LICENSE DESTINATION include/pybind11 COMPONENT ngsolve_devel)
add_subdirectory(netgen)
