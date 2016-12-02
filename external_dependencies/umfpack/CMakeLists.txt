cmake_minimum_required(VERSION 2.8)
project(umfpack)
include(ExternalProject)

message("LAPACK ${LAPACK_LIBRARIES}")
ExternalProject_Add(
  suitesparse
  PREFIX ${CMAKE_CURRENT_BINARY_DIR}/umfpack
  GIT_REPOSITORY https://github.com/jlblancoc/suitesparse-metis-for-windows.git
  CMAKE_ARGS -DSUITESPARSE_USE_CUSTOM_BLAS_LAPACK_LIBS=ON -DSHARED=ON -DBUILD_METIS=OFF -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR} -DSUITESPARSE_CUSTOM_LAPACK_LIB=${LAPACK_LIBRARIES} -DSUITESPARSE_CUSTOM_BLAS_LIB=${LAPACK_LIBRARIES}
  )

