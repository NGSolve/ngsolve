include(CheckCXXSourceCompiles)
include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})

if(EMSCRIPTEN)
  add_custom_command(OUTPUT kernel_generator
    COMMAND g++ -std=c++17 ${CMAKE_CURRENT_SOURCE_DIR}/generate_mat_kernels.cpp -o ${CMAKE_CURRENT_BINARY_DIR}/kernel_generator
    )
  set(kernel_command ${CMAKE_CURRENT_BINARY_DIR}/kernel_generator)
else(EMSCRIPTEN)
  add_executable(kernel_generator generate_mat_kernels.cpp)
  set(kernel_command kernel_generator)
endif()

add_custom_command(OUTPUT  matkernel_1.hpp matkernel_2.hpp matkernel_4.hpp matkernel_8.hpp
  COMMAND ${kernel_command} matkernel_1.hpp 1
  COMMAND ${kernel_command} matkernel_2.hpp 2
  COMMAND ${kernel_command} matkernel_4.hpp 4
  COMMAND ${kernel_command} matkernel_8.hpp 8
  DEPENDS kernel_generator
)

add_custom_target(kernel_generated DEPENDS matkernel_1.hpp matkernel_2.hpp matkernel_4.hpp matkernel_8.hpp)

add_library(ngbla ${NGS_LIB_TYPE}
        bandmatrix.cpp triangular.cpp calcinverse.cpp cholesky.cpp
        LUdecomposition.cpp householder.cpp svd.cpp
        eigensystem.cpp LapackGEP.cpp
        python_bla.cpp ngblas.cpp ngblasc.cpp timing_ngblas.cpp
        )

add_dependencies(ngbla kernel_generated)

target_include_directories(ngbla PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
target_compile_definitions(ngbla PRIVATE ${NGSOLVE_COMPILE_DEFINITIONS_PRIVATE})

target_link_libraries(ngbla PUBLIC ngstd ${MPI_CXX_LIBRARIES} PRIVATE "$<BUILD_INTERFACE:netgen_python>")
target_link_libraries(ngbla ${LAPACK_CMAKE_LINK_INTERFACE} "$<BUILD_INTERFACE:ngs_lapack>")

install( TARGETS ngbla ${ngs_install_dir} )

install( FILES
        bandmatrix.hpp cholesky.hpp matrix.hpp ng_lapack.hpp 
        vector.hpp bla.hpp expr.hpp symmetricmatrix.hpp clapack.h     
        tensor.hpp triangular.hpp householder.hpp eigensystem.hpp calcinverse.hpp
        ngblas.hpp simd_complex.hpp complex_wrapper.hpp
        DESTINATION ${NGSOLVE_INSTALL_DIR_INCLUDE}
        COMPONENT ngsolve_devel
       )

