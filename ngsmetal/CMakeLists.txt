include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})

# metal-cpp (header-only) is downloaded at configure time.
set(METAL_CPP_VERSION 27c4382b7151d55a51692cdcb27aaa98752240de CACHE STRING "metal-cpp git commit/tag to download")

include(FetchContent)
FetchContent_Declare(
  metal_cpp
  URL "https://github.com/apple/metal-cpp/archive/${METAL_CPP_VERSION}.tar.gz"
  URL_HASH SHA256=01b6e76f10c0f5b6eb44be22a3d276cd9922da62a51d4d345b70c5177d720231
  DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../external_dependencies
  DOWNLOAD_EXTRACT_TIMESTAMP FALSE
)
FetchContent_MakeAvailable(metal_cpp)

if (NOT COMMAND pybind11_add_module)
  add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../external_dependencies/netgen/external_dependencies/pybind11 pybind11_build)
endif()

pybind11_add_module(libngsmetal python_ngsmetal.cpp metal_device.cpp )

target_include_directories(libngsmetal PRIVATE ${NETGEN_PYTHON_INCLUDE_DIRS} ${metal_cpp_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/../ngstd ${CMAKE_CURRENT_SOURCE_DIR}/../linalg)

target_link_libraries (libngsmetal PUBLIC ngcomp ngfem ngla ngbla ngstd ${MPI_CXX_LIBRARIES} PRIVATE "$<BUILD_INTERFACE:netgen_python>" ${HYPRE_LIBRARIES})
target_link_libraries(libngsmetal ${LAPACK_CMAKE_LINK_INTERFACE} "$<BUILD_INTERFACE:ngs_lapack>")
install( TARGETS libngsmetal ${ngs_install_dir} )

install( FILES
        gpu_btdtb.hpp metal_device.hpp
        DESTINATION ${NGSOLVE_INSTALL_DIR_INCLUDE}
        COMPONENT ngsolve_devel
       )

