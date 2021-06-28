if(NOT UNIX)
  message(FATAL_ERROR "Building HYPRE as dependency is not supported on this platform. Please configure with USE_HYPRE=OFF or set HYPRE_DIR=path_to_your_hypre_installation")
endif(NOT UNIX)

set(HYPRE_SRC_DIR ${CMAKE_CURRENT_BINARY_DIR}/dependencies/src/project_parmetis/src)
set(HYPRE_DIR ${CMAKE_CURRENT_BINARY_DIR}/dependencies/src/project_hypre/src/hypre)

ExternalProject_Add(project_hypre
  PREFIX ${CMAKE_CURRENT_BINARY_DIR}/dependencies
  URL "https://github.com/hypre-space/hypre/archive/refs/tags/v2.11.1.tar.gz"
  URL_MD5 28f3928b062c79c2eaf54c0978efcbfb
  CONFIGURE_COMMAND cmake src
  -DCMAKE_POSITION_INDEPENDENT_CODE=ON
  -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
  -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
  BUILD_IN_SOURCE 1
  DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external_dependencies
  PATCH_COMMAND ""
  UPDATE_COMMAND "" # Disable update
  )

set_vars( NGSOLVE_CMAKE_ARGS HYPRE_DIR )

list(APPEND DEPENDENCIES project_hypre)
