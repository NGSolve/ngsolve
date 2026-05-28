if(NOT UNIX)
  message(FATAL_ERROR "Building HYPRE as dependency is not supported on this platform. Please configure with USE_HYPRE=OFF or set HYPRE_DIR=path_to_your_hypre_installation")
endif(NOT UNIX)

set(HYPRE_DIR ${CMAKE_CURRENT_BINARY_DIR}/dependencies/hypre)

set(HYPRE_CMAKE_ARGS
  ${SUBPROJECT_CMAKE_ARGS}
  -DCMAKE_INSTALL_PREFIX=${HYPRE_DIR}
  -DHYPRE_WITH_MPI=${USE_MPI}
  -DHYPRE_WITH_OPENMP=ON
)

if(USE_MPI)
  list(APPEND HYPRE_CMAKE_ARGS
    -DMPI_C_COMPILER=${MPI_C_COMPILER}
    -DMPI_CXX_COMPILER=${MPI_CXX_COMPILER}
  )
endif(USE_MPI)

ExternalProject_Add(project_hypre
  ${SUBPROJECT_ARGS}
  URL "https://github.com/hypre-space/hypre/archive/refs/tags/v2.33.0.tar.gz"
  URL_MD5 d4990384b7b1d8b0357fc34d91530d49
  SOURCE_SUBDIR src
  CMAKE_ARGS ${HYPRE_CMAKE_ARGS}
  DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external_dependencies
  PATCH_COMMAND ""
  UPDATE_COMMAND ""
  )

set_vars( NGSOLVE_CMAKE_ARGS HYPRE_DIR )

list(APPEND DEPENDENCIES project_hypre)
