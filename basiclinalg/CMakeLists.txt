include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(kernel_generator generate_mat_kernels.cpp)

target_compile_definitions(kernel_generator PUBLIC ${NGSOLVE_COMPILE_DEFINITIONS})
target_compile_definitions(kernel_generator PRIVATE ${NGSOLVE_COMPILE_DEFINITIONS_PRIVATE})
target_compile_options(kernel_generator PUBLIC ${NGSOLVE_COMPILE_OPTIONS} ${ngcore_compile_options})

add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/matkernel.hpp
  COMMAND kernel_generator
  DEPENDS kernel_generator
  )

add_custom_target(kernel_generated DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/matkernel.hpp)

add_library(ngbla ${NGS_LIB_TYPE}
        bandmatrix.cpp triangular.cpp calcinverse.cpp cholesky.cpp
        LUdecomposition.cpp householder.cpp
        eigensystem.cpp LapackGEP.cpp
        python_bla.cpp avector.cpp ngblas.cpp
        )

add_dependencies(ngbla kernel_generated)

target_include_directories(ngbla PRIVATE ${CMAKE_CURRENT_BINARY_DIR} ${NETGEN_PYTHON_INCLUDE_DIRS})
target_compile_definitions(ngbla PRIVATE ${NGSOLVE_COMPILE_DEFINITIONS_PRIVATE})

if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  target_compile_options(ngbla PUBLIC -Wno-undefined-var-template)
endif()

target_link_libraries(ngbla PUBLIC ngstd ${MPI_CXX_LIBRARIES} PRIVATE netgen_python)
target_link_libraries(ngbla ${LAPACK_CMAKE_LINK_INTERFACE} ${LAPACK_LIBRARIES})

install( TARGETS ngbla ${ngs_install_dir} )

install( FILES
        bandmatrix.hpp cholesky.hpp matrix.hpp ng_lapack.hpp 
        vector.hpp bla.hpp expr.hpp symmetricmatrix.hpp arch.hpp clapack.h     
        tensor.hpp triangular.hpp householder.hpp cuda_bla.hpp avector.hpp ngblas.hpp
        DESTINATION ${NGSOLVE_INSTALL_DIR_INCLUDE}
        COMPONENT ngsolve_devel
       )

