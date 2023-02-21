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

find_package(CUDAToolkit)
enable_language(CUDA)

add_compile_definitions(NGSolve PRIVATE CUDA)

include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})

add_library(${lib_name} ${NGS_LIB_TYPE}
    cuda_linalg.cpp unifiedvector.cpp
    cuda_ngstd.cpp cuda_ngbla.hpp cuda_applyIR.cpp
    linalg_kernels.cu 
)

target_link_libraries(${lib_name} PUBLIC CUDA::cusparse CUDA::cublas CUDA::cudart ngsolve)
set_target_properties(${lib_name} PROPERTIES POSITION_INDEPENDENT_CODE ON )
if(CMAKE_PROJECT_NAME STREQUAL "NGSolve")
  install( TARGETS ${lib_name} ${ngs_install_dir} )
endif()

if(NETGEN_USE_PYTHON)
    find_package(Python3 COMPONENTS Interpreter Development)
    Python3_add_library(ngscuda python_ngscuda.cpp)
    target_link_libraries(ngscuda PUBLIC ${lib_name})
    set_target_properties(ngscuda PROPERTIES INSTALL_RPATH "${NETGEN_RPATH_TOKEN}/../${NETGEN_PYTHON_RPATH}")
    install(TARGETS ngscuda DESTINATION ${NGSOLVE_INSTALL_DIR_PYTHON}/ngsolve COMPONENT ngsolve)
endif(NETGEN_USE_PYTHON)

install( FILES
        cuda_linalg.hpp
        unifiedvector.hpp
        cuda_ngstd.hpp
        DESTINATION ${NGSOLVE_INSTALL_DIR_INCLUDE}
        COMPONENT ngsolve_devel
)

