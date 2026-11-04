include(CheckCXXSourceCompiles)
include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})

add_library(ngbla ${NGS_LIB_TYPE}
        bandmatrix.cpp triangular.cpp calcinverse.cpp cholesky.cpp
        LUdecomposition.cpp householder.cpp svd.cpp
        eigensystem.cpp LapackGEP.cpp
        python_bla.cpp ngblas.cpp ngblas32.cpp ngblasc.cpp timing_ngblas.cpp
        )

target_include_directories(ngbla PUBLIC ${NGSOLVE_INCLUDE_DIRS})
target_compile_definitions(ngbla PUBLIC ${NGSOLVE_COMPILE_DEFINITIONS})
target_compile_definitions(ngbla PRIVATE ${NGSOLVE_COMPILE_DEFINITIONS_PRIVATE})
target_compile_options(ngbla PUBLIC ${NGSOLVE_COMPILE_OPTIONS})

target_link_libraries(ngbla PUBLIC netgen_core_lib ${MPI_CXX_LIBRARIES} PRIVATE "$<BUILD_INTERFACE:netgen_python>")
target_link_libraries(ngbla ${LAPACK_CMAKE_LINK_INTERFACE} "$<BUILD_INTERFACE:ngs_lapack>")

install( TARGETS ngbla ${ngs_install_dir} )

install( FILES
        bandmatrix.hpp cholesky.hpp matrix.hpp ng_lapack.hpp 
        vector.hpp bla.hpp expr.hpp symmetricmatrix.hpp clapack.h     
        tensor.hpp triangular.hpp householder.hpp eigensystem.hpp calcinverse.hpp
        ngblas.hpp simd_complex.hpp complex_wrapper.hpp ectypes.hpp
        DESTINATION ${NGSOLVE_INSTALL_DIR_INCLUDE}
        COMPONENT ngsolve_devel
       )
