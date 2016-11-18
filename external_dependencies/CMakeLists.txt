install(DIRECTORY pybind11/include/pybind11 DESTINATION include COMPONENT ngsolve_devel)
install(FILES pybind11/LICENSE DESTINATION include/pybind11 COMPONENT ngsolve_devel)
add_subdirectory(netgen)
