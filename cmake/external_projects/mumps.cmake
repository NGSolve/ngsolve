if(NOT UNIX)
  message(FATAL_ERROR "Building MUMPS as dependency is not supported on this platform. Please configure with USE_MUMPS=OFF or set MUMPS_DIR=path_to_your_mumps_installation")
endif(NOT UNIX)

enable_language(Fortran)
find_package(MPI REQUIRED)

if(NOT PARMETIS_DIR)
  find_package(PARMETIS REQUIRED)
  if (NOT PARMETIS_DIR)
    message(FATAL_ERROR "Do not have a PARMETIS_DIR for MUMPS!")
  endif(NOT PARMETIS_DIR)
endif(NOT PARMETIS_DIR)

set(MUMPS_SRC_DIR ${CMAKE_CURRENT_BINARY_DIR}/dependencies/src/project_mumps)
set(MUMPS_DIR ${MUMPS_SRC_DIR})

configure_file(${CMAKE_CURRENT_LIST_DIR}/mumps.inc ${CMAKE_CURRENT_BINARY_DIR}/dependencies/Makefile_mumps.inc)

ExternalProject_Add(project_mumps
  DEPENDS project_parmetis
  URL 
  "https://distfiles.macports.org/mumps/MUMPS_5.2.1.tar.gz"
  "https://mirrors.cloud.tencent.com/macports/distfiles/mumps/MUMPS_5.2.1.tar.gz"
  URL_MD5 a4d43b459dc46db984503fbd8526fa69
  ${SUBPROJECT_ARGS}
  DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external_dependencies
  BUILD_IN_SOURCE 1
  UPDATE_COMMAND cp ${CMAKE_CURRENT_BINARY_DIR}/dependencies/Makefile_mumps.inc ${MUMPS_SRC_DIR}/Makefile.inc
  CONFIGURE_COMMAND ""
  BUILD_COMMAND make alllib -j1
  INSTALL_COMMAND ""
  )

set_vars( NGSOLVE_CMAKE_ARGS MUMPS_DIR )

list(APPEND DEPENDENCIES project_mumps)
