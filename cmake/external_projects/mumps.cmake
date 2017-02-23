set(MUMPS_SRC_DIR ${CMAKE_CURRENT_BINARY_DIR}/dependencies/src/project_mumps)
set(MUMPS_DIR ${MUMPS_SRC_DIR})

find_package(MPI REQUIRED)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cmake/external_projects/mumps.inc ${CMAKE_CURRENT_BINARY_DIR}/Makefile_mumps.inc)

ExternalProject_Add(project_mumps
  DEPENDS project_parmetis
  PREFIX ${CMAKE_CURRENT_BINARY_DIR}/dependencies
#   URL "http://mumps.enseeiht.fr/MUMPS_5.0.2.tar.gz"
#   URL_MD5 591bcb2c205dcb0283872608cdf04927
  URL "http://mumps.enseeiht.fr/MUMPS_5.0.1.tar.gz"
  URL_MD5 b477573fdcc87babe861f62316833db0
  DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external_dependencies
  BUILD_IN_SOURCE 1
  UPDATE_COMMAND cp ${CMAKE_CURRENT_BINARY_DIR}/Makefile_mumps.inc ${MUMPS_SRC_DIR}/Makefile.inc
  CONFIGURE_COMMAND ""
  BUILD_COMMAND make alllib -j1
  INSTALL_COMMAND ""
  )

set_vars( NGSOLVE_CMAKE_ARGS MUMPS_DIR )

list(APPEND DEPENDENCIES project_mumps)
