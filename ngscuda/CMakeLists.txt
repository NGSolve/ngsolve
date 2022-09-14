find_package(CUDAToolkit)
enable_language(CUDA)

add_compile_definitions(NGSolve PRIVATE CUDA)

include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})

add_library(ngscudalib ${NGS_LIB_TYPE}
    cuda_linalg.cpp
		cuda_ngstd.cpp
    linalg_kernels.cu
)

target_link_libraries(ngscudalib PUBLIC CUDA::cusparse CUDA::cublas CUDA::cudart ngsolve)

# TODO: remove. just for testing
message(STATUS "NETGEN_USE_PYTHON VALUE: ${NETGEN_USE_PYTHON}")
message(STATUS "USE_PYTHON VALUE: ${USE_PYTHON}")

if(USE_PYTHON)
    find_package(Python3 COMPONENTS Interpreter Development)
    Python3_add_library(ngscuda python_ngscuda.cpp)
    target_link_libraries(ngscuda PUBLIC ngscudalib)
endif(USE_PYTHON)

install( TARGETS ngscudalib ${ngs_install_dir} )
install(TARGETS ngscuda DESTINATION ${NGSOLVE_INSTALL_DIR_PYTHON}/ngsolve COMPONENT ngsolve)

install( FILES
        cuda_linalg.hpp
				cuda_ngstd.hpp
        DESTINATION ${NGSOLVE_INSTALL_DIR_INCLUDE}
        COMPONENT ngsolve_devel
)

