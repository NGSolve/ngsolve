include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})

add_library(ngla ${NGS_LIB_TYPE}
        basematrix.cpp basevector.cpp multivector.cpp
        blockjacobi.cpp cg.cpp chebyshev.cpp commutingAMG.cpp eigen.cpp	     
        jacobi.cpp order.cpp pardisoinverse.cpp sparsecholesky.cpp	     
        sparsematrix.cpp sparsematrix_dyn.cpp special_matrix.cpp superluinverse.cpp		     
        mumpsinverse.cpp elementbyelement.cpp arnoldi.cpp paralleldofs.cpp diagonalmatrix.cpp
        python_linalg.cpp umfpackinverse.cpp
        ../parallel/parallelvvector.cpp ../parallel/parallel_matrices.cpp 
)

target_include_directories(ngla PRIVATE ${UMFPACK_INCLUDE_DIR} ${NETGEN_PYTHON_INCLUDE_DIRS})
target_compile_definitions(ngla PRIVATE ${NGSOLVE_COMPILE_DEFINITIONS_PRIVATE})

if(UMFPACK_STATIC)
    # work-around for private linking of interface libraries
    # see https://gitlab.kitware.com/cmake/cmake/-/issues/15415
    target_link_libraries(ngla PRIVATE $<BUILD_INTERFACE:${UMFPACK_LIBRARIES}>)
else(UMFPACK_STATIC)
    target_link_libraries(ngla PRIVATE ${UMFPACK_LIBRARIES})
endif()

target_link_libraries(ngla ${LAPACK_CMAKE_LINK_INTERFACE} "$<BUILD_INTERFACE:ngs_lapack>")

target_link_libraries(ngla PUBLIC ngbla ngstd PRIVATE "$<BUILD_INTERFACE:netgen_python>"
    ${MUMPS_LIBRARIES} ${MPI_Fortran_LIBRARIES} ${MPI_CXX_LIBRARIES})
install( TARGETS ngla ${ngs_install_dir} )

if(USE_MUMPS)
    # Create dummy fortran lib
    # This forces cmake to link fortran runtime libraries to ngla
    file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/dummy.f "")
    add_library(dummy_fortran ${CMAKE_CURRENT_BINARY_DIR}/dummy.f)
    target_link_libraries(ngla PUBLIC dummy_fortran )
    install( TARGETS dummy_fortran ${ngs_install_dir} )
endif()

install( FILES
        basematrix.hpp basevector.hpp blockjacobi.hpp cg.hpp multivector.hpp
        chebyshev.hpp commutingAMG.hpp eigen.hpp jacobi.hpp la.hpp order.hpp   
        pardisoinverse.hpp sparsecholesky.hpp sparsematrix.hpp
        sparsematrix_impl.hpp sparsematrix_dyn.hpp diagonalmatrix.hpp
        special_matrix.hpp superluinverse.hpp mumpsinverse.hpp
        umfpackinverse.hpp vvector.hpp python_linalg.hpp
        elementbyelement.hpp arnoldi.hpp paralleldofs.hpp
        DESTINATION ${NGSOLVE_INSTALL_DIR_INCLUDE}
        COMPONENT ngsolve_devel
       )

