if(NOT UNIX)
  message(FATAL_ERROR "Building MUMPS as dependency is not supported on this platform. Please configure with USE_MUMPS=OFF or set MUMPS_DIR=path_to_your_mumps_installation")
endif(NOT UNIX)

if(NOT PARMETIS_DIR)
  include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/external_projects/parmetis.cmake)
endif(NOT PARMETIS_DIR)

set(MUMPS_SRC_DIR ${CMAKE_CURRENT_BINARY_DIR}/dependencies/src/project_mumps)
set(MUMPS_DIR ${MUMPS_SRC_DIR})

enable_language(Fortran)
find_package(MPI REQUIRED)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cmake/external_projects/mumps.inc ${CMAKE_CURRENT_BINARY_DIR}/dependencies/Makefile_mumps.inc)

ExternalProject_Add(project_mumps
  DEPENDS project_parmetis
  PREFIX ${CMAKE_CURRENT_BINARY_DIR}/dependencies
#   URL "http://mumps.enseeiht.fr/MUMPS_5.0.2.tar.gz"
#   URL_MD5 591bcb2c205dcb0283872608cdf04927
  URL "http://mumps.enseeiht.fr/MUMPS_5.0.1.tar.gz"
  URL_MD5 b477573fdcc87babe861f62316833db0
  DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external_dependencies
  BUILD_IN_SOURCE 1
  UPDATE_COMMAND cp ${CMAKE_CURRENT_BINARY_DIR}/dependencies/Makefile_mumps.inc ${MUMPS_SRC_DIR}/Makefile.inc
  CONFIGURE_COMMAND ""
  BUILD_COMMAND make alllib -j1
  INSTALL_COMMAND ""
  )

set_vars( NGSOLVE_CMAKE_ARGS MUMPS_DIR )

list(APPEND DEPENDENCIES project_mumps)
