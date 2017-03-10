include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})

if(NOT WIN32)
    set(ngla_object_libs $<TARGET_OBJECTS:parallel> )
endif(NOT WIN32)
add_library(ngla ${NGS_LIB_TYPE}
        ${ngla_object_libs}
        linalg_kernels.cu basematrix.cpp basevector.cpp 
        blockjacobi.cpp cg.cpp chebyshev.cpp commutingAMG.cpp eigen.cpp	     
        jacobi.cpp order.cpp pardisoinverse.cpp sparsecholesky.cpp	     
        sparsematrix.cpp special_matrix.cpp superluinverse.cpp		     
        mumpsinverse.cpp elementbyelement.cpp arnoldi.cpp paralleldofs.cpp   
        python_linalg.cpp umfpackinverse.cpp
        )

if(NOT WIN32)
  target_link_libraries(ngla ngbla ngstd ${PYTHON_LIBS} ${LAPACK_LIBRARIES} ${MUMPS_LIBRARIES} ${UMFPACK_LIBRARIES} ${SCALAPACK_LIBRARY} ${MPI_Fortran_LIBRARIES} ${MPI_CXX_LIBRARIES})
    install( TARGETS ngla ${ngs_install_dir} )
endif(NOT WIN32)

if(USE_CUDA)
    cuda_add_library (linalg_kernels linalg_kernels.cu cuda_linalg.cpp )
    target_link_libraries(ngla linalg_kernels ${CUDA_LIBRARIES} ${CUDA_CUBLAS_LIBRARIES} ${CUDA_cusparse_LIBRARY})
endif(USE_CUDA)

install( FILES
        basematrix.hpp basevector.hpp blockjacobi.hpp cg.hpp 
        chebyshev.hpp commutingAMG.hpp eigen.hpp jacobi.hpp la.hpp order.hpp   
        pardisoinverse.hpp sparsecholesky.hpp sparsematrix.hpp sparsematrix_spec.hpp
        special_matrix.hpp superluinverse.hpp mumpsinverse.hpp
        umfpackinverse.hpp vvector.hpp     
        elementbyelement.hpp arnoldi.hpp paralleldofs.hpp cuda_linalg.hpp
        DESTINATION include
        COMPONENT ngsolve_devel
       )

