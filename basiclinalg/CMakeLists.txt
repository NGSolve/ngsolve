add_library(ngbla ${NGS_LIB_TYPE}
        bandmatrix.cpp calcinverse.cpp cholesky.cpp 
        eigensystem.cpp vecmat.cpp LapackGEP.cpp LapackInterface.hpp	  
        python_bla.cpp
        )

if(NOT WIN32)
    TARGET_LINK_LIBRARIES(ngbla ngstd ${MPI_CXX_LIBRARIES} ${PYTHON_LIBS} ${LAPACK_LIBRARIES})
    INSTALL( TARGETS ngbla ${ngs_install_dir})
endif(NOT WIN32)

install( FILES
        bandmatrix.hpp cholesky.hpp matrix.hpp ng_lapack.hpp 
        vector.hpp bla.hpp expr.hpp symmetricmatrix.hpp arch.hpp clapack.h     
        tensor.hpp cuda_bla.hpp
        DESTINATION include
       )

