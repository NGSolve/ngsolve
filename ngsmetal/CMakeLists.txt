include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})


if (NOT COMMAND pybind11_add_module)
  add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../external_dependencies/netgen/external_dependencies/pybind11 pybind11_build)
endif()

pybind11_add_module(libngsmetal python_ngsmetal.cpp ngsmetal.cpp metal_btdtb.cpp metal_benchmarks.cpp metal_device.cpp )

target_include_directories(libngsmetal PRIVATE ${NETGEN_PYTHON_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}/metal-cpp ${CMAKE_CURRENT_SOURCE_DIR}/../ngstd ${CMAKE_CURRENT_SOURCE_DIR}/../linalg)

target_link_libraries (libngsmetal PUBLIC ngcomp ngfem ngla ngbla ngstd ${MPI_CXX_LIBRARIES} PRIVATE "$<BUILD_INTERFACE:netgen_python>" ${HYPRE_LIBRARIES})
target_link_libraries(libngsmetal ${LAPACK_CMAKE_LINK_INTERFACE} "$<BUILD_INTERFACE:ngs_lapack>")
install( TARGETS libngsmetal ${ngs_install_dir} )

install( FILES
        ngsmetal.hpp metal_btdtb.hpp metal_device.hpp
        DESTINATION ${NGSOLVE_INSTALL_DIR_INCLUDE}
        COMPONENT ngsolve_devel
       )

