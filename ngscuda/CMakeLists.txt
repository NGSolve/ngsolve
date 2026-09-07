if(NOT CMAKE_PROJECT_NAME STREQUAL "NGSolve")
  cmake_minimum_required(VERSION 3.18)
  project(ngscuda)
  find_package(NGSolve)
  if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX "${NGSOLVE_INSTALL_DIR}" CACHE PATH "Install directory" FORCE)
  endif()
  set(lib_name ngscudalib_local)
else()
  set(lib_name ngscudalib)
endif()

# host compiler only: the kernels are compiled at run time by nvrtc
find_package(CUDAToolkit REQUIRED)

include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})

add_library(${lib_name} SHARED
    cuda_linalg.cpp cuda_ngstd.cpp cuda_ngbla.hpp
    cuda_profiler.cpp cuda_device.cpp
)

target_link_libraries(${lib_name} PUBLIC CUDA::cudart CUDA::nvrtc CUDA::cuda_driver ngsolve $<BUILD_INTERFACE:netgen_python>)
set_target_properties(${lib_name} PROPERTIES POSITION_INDEPENDENT_CODE ON)
if(CMAKE_PROJECT_NAME STREQUAL "NGSolve")
  install( TARGETS ${lib_name} ${ngs_install_dir} )
endif()

if(NETGEN_USE_PYTHON)
    find_package(Python3 COMPONENTS Interpreter Development)
    Python3_add_library(_ngscuda python_ngscuda.cpp)
    target_link_libraries(_ngscuda PUBLIC ${lib_name})
    set_target_properties(_ngscuda PROPERTIES INSTALL_RPATH "${NETGEN_RPATH_TOKEN}/../${NETGEN_PYTHON_RPATH}")
    install(TARGETS _ngscuda DESTINATION ${NGSOLVE_INSTALL_DIR_PYTHON}/ngsolve COMPONENT ngsolve)
endif(NETGEN_USE_PYTHON)

install( FILES
        cuda_linalg.hpp
        cuda_profiler.hpp
        cuda_core.hpp
        cuda_ngstd.hpp
        cuda_ngbla.hpp
        cuda_device.hpp
        DESTINATION ${NGSOLVE_INSTALL_DIR_INCLUDE}
        COMPONENT ngsolve_devel
)

