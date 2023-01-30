find_package(CUDAToolkit)
enable_language(CUDA)

add_compile_definitions(NGSolve PRIVATE CUDA)

include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})

add_library(ngscudalib ${NGS_LIB_TYPE}
    cuda_linalg.cpp unifiedvector.cpp
    cuda_ngstd.cpp cude_applyIR.cpp
    linalg_kernels.cu
)

target_link_libraries(ngscudalib PUBLIC CUDA::cusparse CUDA::cublas CUDA::cudart ngsolve)
install( TARGETS ngscudalib ${ngs_install_dir} )

if(NETGEN_USE_PYTHON)
    find_package(Python3 COMPONENTS Interpreter Development)
    Python3_add_library(ngscuda python_ngscuda.cpp)
    target_link_libraries(ngscuda PUBLIC ngscudalib)
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

