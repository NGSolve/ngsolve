include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})

add_library(ngbla ${NGS_LIB_TYPE}
        bandmatrix.cpp calcinverse.cpp cholesky.cpp 
        eigensystem.cpp vecmat.cpp LapackGEP.cpp
        python_bla.cpp avector.cpp
        )

if(NOT WIN32)
    TARGET_LINK_LIBRARIES(ngbla ngstd ${MPI_CXX_LIBRARIES} ${NETGEN_PYTHON_LIBRARIES} ${LAPACK_LIBRARIES})
    INSTALL( TARGETS ngbla ${ngs_install_dir} )
endif(NOT WIN32)

install( FILES
        bandmatrix.hpp cholesky.hpp matrix.hpp ng_lapack.hpp 
        vector.hpp bla.hpp expr.hpp symmetricmatrix.hpp arch.hpp clapack.h     
        tensor.hpp cuda_bla.hpp avector.hpp
        DESTINATION include
        COMPONENT ngsolve_devel
       )

