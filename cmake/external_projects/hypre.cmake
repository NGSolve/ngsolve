if(NOT UNIX)
  message(FATAL_ERROR "Building HYPRE as dependency is not supported on this platform. Please configure with USE_HYPRE=OFF or set HYPRE_DIR=path_to_your_hypre_installation")
endif(NOT UNIX)

set(HYPRE_SRC_DIR ${CMAKE_CURRENT_BINARY_DIR}/dependencies/src/project_parmetis/src)
set(HYPRE_DIR ${CMAKE_CURRENT_BINARY_DIR}/dependencies/src/project_hypre/src/hypre)

ExternalProject_Add(project_hypre
  PREFIX ${CMAKE_CURRENT_BINARY_DIR}/dependencies
  URL "http://computation.llnl.gov/projects/hypre-scalable-linear-solvers-multigrid-methods/download/hypre-2.11.1.tar.gz"
  URL_MD5 3f02ef8fd679239a6723f60b7f796519
  CONFIGURE_COMMAND cmake src -DCMAKE_C_FLAGS=-fPIC
  BUILD_IN_SOURCE 1
  DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external_dependencies
  PATCH_COMMAND ""
  UPDATE_COMMAND "" # Disable update
  )

set_vars( NGSOLVE_CMAKE_ARGS HYPRE_DIR )

list(APPEND DEPENDENCIES project_hypre)
