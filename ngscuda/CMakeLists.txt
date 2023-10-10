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

add_compile_definitions(CUDA)

include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})

add_library(${lib_name} SHARED
    cuda_linalg.cpp unifiedvector.cpp
    cuda_ngstd.cpp cuda_ngbla.hpp cuda_applyIR.cpp
    linalg_kernels.cu dev_sparsecholesky.cpp dev_blockjacobi.cpp cuda_profiler.cu
)

set_source_files_properties(cuda_linalg.cpp unifiedvector.cpp dev_sparsecholesky.cpp dev_blockjacobi.cpp PROPERTIES LANGUAGE CUDA)

target_compile_options(${lib_name} PUBLIC $<$<COMPILE_LANGUAGE:CUDA>:--expt-relaxed-constexpr --extended-lambda >)

target_link_libraries(${lib_name} PUBLIC CUDA::cusparse CUDA::cublas CUDA::cudart ngsolve netgen_python)
set_target_properties(${lib_name} PROPERTIES POSITION_INDEPENDENT_CODE ON CUDA_SEPARABLE_COMPILATION ON)
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

set( ngs_nvcc_header "\
${ngscxx_set_script_dir}
if [ -f ${CMAKE_CUDA_COMPILER} ]
  then NVCC=${CMAKE_CUDA_COMPILER}
  else NVCC=nvcc
fi
")

set(nvcc_flags "")
foreach( opt ${NGSOLVE_COMPILE_OPTIONS} ${ngcore_compile_options} "-fPIC")
    set(nvcc_flags "${nvcc_flags} -Xcompiler ${opt}")
    message("add opt ", ${opt})
endforeach()

file(GENERATE OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/ngs_nvcc
              CONTENT "${ngs_nvcc_header}\n$NVCC ${ngscxx_define_flags} -arch=native ${nvcc_flags} ${ngscxx_includes} $*\n"
              CONDITION $<COMPILE_LANGUAGE:CXX>
              )
file(GENERATE OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/ngs_nvlink CONTENT "${ngs_nvcc_header}\n$NGSCXX ${ngsld_flags} $*\n")
install (PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/ngs_nvcc DESTINATION ${NGSOLVE_INSTALL_DIR_BIN} COMPONENT ngsolve_devel )
install (PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/ngs_nvlink DESTINATION ${NGSOLVE_INSTALL_DIR_BIN} COMPONENT ngsolve_devel )

install( FILES
        cuda_linalg.hpp
        unifiedvector.hpp
        cuda_ngstd.hpp
        DESTINATION ${NGSOLVE_INSTALL_DIR_INCLUDE}
        COMPONENT ngsolve_devel
)

