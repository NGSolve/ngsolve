include (ExternalProject)

set_property (DIRECTORY PROPERTY EP_BASE dependencies)

set (DEPENDENCIES)
set (LAPACK_PROJECTS)
set (NGSOLVE_CMAKE_ARGS)
set (NETGEN_CMAKE_ARGS)

set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH}" "${PROJECT_SOURCE_DIR}/cmake_modules")

macro(set_vars OUTPUT_VAR )
  foreach(varname ${ARGN})
    if(NOT "${${varname}}" STREQUAL "")
      string(REPLACE ";" "$<SEMICOLON>" varvalue "${${varname}}" )
      list(APPEND ${OUTPUT_VAR} -D${varname}=${varvalue})
    endif()
  endforeach()
endmacro()
#######################################################################
if(WIN32)
  if(NOT CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
    string(REGEX REPLACE "/W[0-4]" "/W0" CMAKE_CXX_FLAGS_NEW ${CMAKE_CXX_FLAGS})
    set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS_NEW} CACHE STRING "compile flags" FORCE)
    string(REGEX REPLACE "/W[0-4]" "/W0" CMAKE_CXX_FLAGS_NEW ${CMAKE_CXX_FLAGS_RELEASE})
    set(CMAKE_CXX_FLAGS_RELEASE ${CMAKE_CXX_FLAGS_NEW} CACHE STRING "compile flags" FORCE)

    string(REGEX REPLACE "/W[0-4]" "/W0" CMAKE_SHARED_LINKER_FLAGS_NEW ${CMAKE_SHARED_LINKER_FLAGS})
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS_NEW} /IGNORE:4217,4049" CACHE STRING "compile flags" FORCE)
    string(REGEX REPLACE "/W[0-4]" "/W0" CMAKE_EXE_LINKER_FLAGS_NEW ${CMAKE_EXE_LINKER_FLAGS})
    set(CMAKE_EXE_LINKER_FLAGS"${CMAKE_EXE_LINKER_FLAGS_NEW}/IGNORE:4217,4049" CACHE STRING "compile flags" FORCE)

    set_vars(NGSOLVE_CMAKE_ARGS CMAKE_SHARED_LINKER_FLAGS CMAKE_SHARED_LINKER_FLAGS_RELEASE CMAKE_CXX_FLAGS CMAKE_CXX_FLAGS_RELEASE)
  endif(NOT CMAKE_CXX_COMPILER_ID STREQUAL "Intel")

  if(${CMAKE_SIZEOF_VOID_P} MATCHES 4)
    # 32 bit
    set(LAPACK_DOWNLOAD_URL_WIN "http://www.asc.tuwien.ac.at/~mhochsteger/ngsuite/lapack32.zip" CACHE STRING INTERNAL)
  else(${CMAKE_SIZEOF_VOID_P} MATCHES 4)
    # 64 bit
    set(LAPACK_DOWNLOAD_URL_WIN "http://www.asc.tuwien.ac.at/~mhochsteger/ngsuite/lapack64.zip" CACHE STRING INTERNAL)
  endif(${CMAKE_SIZEOF_VOID_P} MATCHES 4)
endif(WIN32)

#######################################################################
option(BUILD_NETGEN ON "Build Netgen from the git submodule as part of NGSolve (recommended)")
#######################################################################
if(NETGEN_SOURCE_DIR)
  message(STATUS "Since NETGEN_SOURCE_DIR is given, assume Netgen is already installed")
  message(STATUS "Looking for NetgenConfig.cmake...")
  find_package(Netgen REQUIRED CONFIG HINTS ${INSTALL_DIR}/share/cmake $ENV{NETGENDIR}/../share/cmake)
  set(INSTALL_DIR ${NETGEN_INSTALL_DIR} CACHE PATH "Install path")
  set(BUILD_NETGEN OFF)
else(NETGEN_SOURCE_DIR)
  message(STATUS "Build Netgen from git submodule")
#   execute_process(COMMAND git submodule update --init --recursive WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
  execute_process(COMMAND cmake -P cmake_modules/check_submodules.cmake WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} )
  add_custom_target(check_submodules_start ALL cmake -P cmake_modules/check_submodules.cmake WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} )
  add_custom_target(check_submodules_stop ALL cmake -P cmake_modules/check_submodules.cmake WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} DEPENDS ngsolve)
  set(BUILD_NETGEN ON)
  set(INSTALL_DIR /opt/netgen CACHE PATH "Install path")
  set(NETGEN_DIR ${INSTALL_DIR})
  set (NETGEN_CMAKE_ARGS)
  set_vars( NETGEN_CMAKE_ARGS
    USE_GUI
    USE_PYTHON
    USE_MPI
    USE_CCACHE
    USE_NATIVE_ARCH
    USE_OCC
    INSTALL_DIR
    INSTALL_DEPENDENCIES 
    INTEL_MIC 
  )
endif(NETGEN_SOURCE_DIR)

set_vars(NGSOLVE_CMAKE_ARGS NETGEN_DIR)

#######################################################################
# find netgen

#######################################################################
set(LAPACK_LIBRARIES CACHE INTERNAL "Lapack libraries")
if(USE_MKL)
    set(MKL_MULTI_THREADED ON)

    set(MKL_STATIC OFF CACHE BOOL "Link static MKL")
    set(MKL_SDL ON CACHE BOOL "Link single dynamic MKL lib")

    set(USE_LAPACK ON)
    find_package(MKL REQUIRED)
    if(USE_MUMPS)
        # include scalapack
        set( LAPACK_LIBRARIES "${MKL_LIBRARIES}")
    else(USE_MUMPS)
        set( LAPACK_LIBRARIES "${MKL_MINIMAL_LIBRARIES}")
    endif(USE_MUMPS)
    set_vars(NGSOLVE_CMAKE_ARGS MKL_LINK_FLAGS MKL_INTERFACE_LAYER MKL_INCLUDE_DIRS)
endif(USE_MKL)

if (USE_LAPACK)
    if(NOT LAPACK_LIBRARIES)
      if(WIN32)
        ExternalProject_Add(win_download_lapack
          PREFIX ${CMAKE_CURRENT_BINARY_DIR}/tcl
          URL ${LAPACK_DOWNLOAD_URL_WIN}
          UPDATE_COMMAND "" # Disable update
          BUILD_IN_SOURCE 1
          CONFIGURE_COMMAND ""
          BUILD_COMMAND ""
          INSTALL_COMMAND ${CMAKE_COMMAND} -E copy_directory . ${INSTALL_DIR}
          )
        set(LAPACK_LIBRARIES ${INSTALL_DIR}/lib/BLAS.lib)
        list(APPEND LAPACK_PROJECTS win_download_lapack)
      else(WIN32)
        find_package(LAPACK)
      endif(WIN32)
    endif()
    set_vars(NGSOLVE_CMAKE_ARGS LAPACK_LIBRARIES)
endif (USE_LAPACK)

#######################################################################

if(USE_UMFPACK)
  ExternalProject_Add(
    suitesparse
    DEPENDS ${LAPACK_PROJECTS}
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/umfpack
    GIT_REPOSITORY https://github.com/jlblancoc/suitesparse-metis-for-windows.git
    CMAKE_ARGS -DSUITESPARSE_USE_CUSTOM_BLAS_LAPACK_LIBS=ON -DSHARED=ON -DBUILD_METIS=OFF -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR} -DSUITESPARSE_CUSTOM_LAPACK_LIB=${LAPACK_LIBRARIES} -DSUITESPARSE_CUSTOM_BLAS_LIB=${LAPACK_LIBRARIES}
    )
  list(APPEND DEPENDENCIES suitesparse)
endif(USE_UMFPACK)

#######################################################################
# propagate cmake variables to NGSolve subproject
set_vars( NGSOLVE_CMAKE_ARGS
  USE_GUI
  USE_PYTHON
  USE_LAPACK
  USE_MPI
  USE_VT
  USE_CUDA
  USE_MKL
  USE_MUMPS
  USE_PARDISO
  USE_UMFPACK
  USE_VTUNE
  USE_NUMA
  USE_CCACHE
  USE_NATIVE_ARCH
  USE_OCC
  INSTALL_DIR
  NETGEN_SOURCE_DIR
  INSTALL_DEPENDENCIES 
  INTEL_MIC 
  )

if(BUILD_NETGEN)
  ExternalProject_Add (netgen_project
    SOURCE_DIR ${PROJECT_SOURCE_DIR}/external_dependencies/netgen
    CMAKE_ARGS ${NETGEN_CMAKE_ARGS}
    BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/netgen
  )
#   list(APPEND DEPENDENCIES netgen_project)
endif(BUILD_NETGEN)

ExternalProject_Add (ngsolve
  DEPENDS ${DEPENDENCIES} ${LAPACK_PROJECTS}
  SOURCE_DIR ${PROJECT_SOURCE_DIR}
  CMAKE_ARGS ${NGSOLVE_CMAKE_ARGS} -DUSE_SUPERBUILD=OFF -DCMAKE_PREFIX_PATH=${NETGEN_DIR}
  INSTALL_COMMAND ""
  BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/ngsolve
  BUILD_COMMAND ${CMAKE_COMMAND} --build ${CMAKE_CURRENT_BINARY_DIR}/ngsolve --config ${CMAKE_BUILD_TYPE}
  STEP_TARGETS build
)
add_dependencies(ngsolve netgen_project)

install(CODE "execute_process(COMMAND cmake --build . --target install WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/ngsolve)")

add_custom_target(test_ngsolve
  ${CMAKE_COMMAND} --build ${CMAKE_CURRENT_BINARY_DIR}/ngsolve
                   --target test
                   --config ${CMAKE_BUILD_TYPE}
                   )
