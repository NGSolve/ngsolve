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
    cuda_profiler.cpp cuda_matrix_free_btdtb.cpp cuda_device.cpp
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

# ngs_nvcc / ngs_nvlink: wrappers for users compiling their own cuda code
# against ngsolve, generated only if a cuda compiler is around
find_program(NGS_NVCC nvcc HINTS ${CUDAToolkit_BIN_DIR})
if(NGS_NVCC)
  set( ngs_nvcc_header "${ngscxx_header}
NVCC=${NGS_NVCC}
")
  set(nvcc_flags "--expt-relaxed-constexpr --extended-lambda -rdc=true")
  foreach( opt ${NGSOLVE_COMPILE_OPTIONS} ${ngcore_compile_options} "-fPIC")
    if(opt MATCHES "-std=")
      set(nvcc_flags "${nvcc_flags} ${opt}")
    else()
      set(nvcc_flags "${nvcc_flags} -Xcompiler ${opt}")
    endif()
  endforeach()
  set(ngs_nvcc_command "$NVCC ${ngscxx_define_flags} -arch=native ${nvcc_flags} ${ngscxx_includes} $*\n")
  if(CCACHE_FOUND)
    set(ngs_nvcc_command "
if [ -f ${CCACHE_FOUND} ]
  then ${CCACHE_FOUND} ${ngs_nvcc_command}
  else ${ngs_nvcc_command}
fi")
  endif()
  set(ngs_nvlink_command "$NVCC -rdc=true -arch=native $*\n")
  file(GENERATE OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/ngs_nvcc CONTENT "${ngs_nvcc_header}\n${ngs_nvcc_command}\n"
       CONDITION $<COMPILE_LANGUAGE:CXX>)
  file(GENERATE OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/ngs_nvlink CONTENT "${ngs_nvcc_header}\n${ngs_nvlink_command}\n")
  install (PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/ngs_nvcc DESTINATION ${NGSOLVE_INSTALL_DIR_BIN} COMPONENT ngsolve_devel )
  install (PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/ngs_nvlink DESTINATION ${NGSOLVE_INSTALL_DIR_BIN} COMPONENT ngsolve_devel )
endif()

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

