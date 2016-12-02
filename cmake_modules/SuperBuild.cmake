include (ExternalProject)

set_property (DIRECTORY PROPERTY EP_BASE Dependencies)

set (DEPENDENCIES)
set (LAPACK_DEPENDENCIES)
set (EXTRA_CMAKE_ARGS)

set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH}" "${PROJECT_SOURCE_DIR}/cmake_modules")

macro(set_vars ...)
  foreach(varname ${ARGV})
    if(NOT "${${varname}}" STREQUAL "")
      string(REPLACE ";" "$<SEMICOLON>" varvalue "${${varname}}" )
      list(APPEND EXTRA_CMAKE_ARGS -D${varname}=${varvalue})
    endif()
  endforeach()
endmacro()

#######################################################################
# find netgen
set(INSTALL_DIR /opt/netgen CACHE PATH "Install path")
if(APPLE)
  set(CMAKE_INSTALL_PREFIX "${INSTALL_DIR}/Netgen.app/Contents/Resources" CACHE INTERNAL "Prefix prepended to install directories" FORCE)
else(APPLE)
  set(CMAKE_INSTALL_PREFIX "${INSTALL_DIR}" CACHE INTERNAL "Prefix prepended to install directories" FORCE) 
endif(APPLE)

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
    set_vars(MKL_LINK_FLAGS MKL_INTERFACE_LAYER MKL_INCLUDE_DIRS)
endif(USE_MKL)

if (USE_LAPACK)
    if(NOT LAPACK_LIBRARIES)
      if(WIN32)
        ExternalProject_Add(win_download_lapack
          PREFIX ${CMAKE_CURRENT_BINARY_DIR}/tcl
          URL "http://www.asc.tuwien.ac.at/~mhochsteger/ngsuite/lapack64.zip"
          UPDATE_COMMAND "" # Disable update
          BUILD_IN_SOURCE 1
          CONFIGURE_COMMAND ""
          BUILD_COMMAND ""
          INSTALL_COMMAND ${CMAKE_COMMAND} -E copy_directory . ${INSTALL_DIR}
          )
        set(LAPACK_LIBRARIES ${INSTALL_DIR}/lib/BLAS.lib)
        list(APPEND LAPACK_DEPENDENCIES win_download_lapack)
      else(WIN32)
        find_package(LAPACK)
      endif(WIN32)
    endif()
    set_vars(LAPACK_LIBRARIES)
endif (USE_LAPACK)

#######################################################################
set(NETGEN_INTERNAL_INCLUDE_DIR NOTFOUND CACHE INTERNAL "Netgen include directories")
set(BUILD_NETGEN OFF CACHE INTERNAL "is set if netgen is built as external project")
#######################################################################
if(NETGEN_SOURCE_DIR)
    find_path(NETGEN_INCLUDE_DIR nginterface_v2_impl.hpp PATHS ${NETGEN_SOURCE_DIR} ${NETGEN_SOURCE_DIR}/libsrc/include)
    if (NOT NETGEN_INCLUDE_DIR)
      message(FATAL_ERROR "Could not find Netgen source files in NETGEN_SOURCE_DIR, which was set to ${NETGEN_SOURCE_DIR}")
    endif()
    set(NETGEN_INTERNAL_INCLUDE_DIR
      ${NETGEN_INCLUDE_DIR}
      ${NETGEN_INCLUDE_DIR}/../general
      ${NETGEN_INCLUDE_DIR}/../visualization
    )
else(NETGEN_SOURCE_DIR)
  message(STATUS "Use Netgen from submodule, updating modules...")
  execute_process(COMMAND git submodule update --init --recursive WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
  add_custom_target(check_submodules_start ALL cmake -P cmake_modules/check_submodules.cmake WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} )
  add_custom_target(check_submodules_stop ALL cmake -P cmake_modules/check_submodules.cmake WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} DEPENDS ngsolve)
  set(BUILD_NETGEN ON)
  set(NETGEN_INTERNAL_INCLUDE_DIR
    ${CMAKE_CURRENT_SOURCE_DIR}/external_dependencies/netgen/libsrc/include
    ${CMAKE_CURRENT_SOURCE_DIR}/external_dependencies/netgen/libsrc/general
    ${CMAKE_CURRENT_SOURCE_DIR}/external_dependencies/netgen/libsrc/visualization
  )
endif(NETGEN_SOURCE_DIR)

set_vars(NETGEN_INTERNAL_INCLUDE_DIR BUILD_NETGEN)

#######################################################################
if (USE_PYTHON)
  find_path(PYBIND_INCLUDE_DIR pybind11/pybind11.h ${NETGEN_INTERNAL_INCLUDE_DIR}/../../external_dependencies/pybind11/include)
    if( NOT PYBIND_INCLUDE_DIR )
      message(FATAL_ERROR "Could NOT find pybind11!")
    endif( NOT PYBIND_INCLUDE_DIR )
    message(STATUS "Found Pybind11: ${PYBIND_INCLUDE_DIR}")
    set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH}" "${PROJECT_SOURCE_DIR}/cmake_modules/python")
    set(PYTHON_VERSION "3" CACHE STRING "Python version (only Python >= 3.0 supported)")
    set(Python_ADDITIONAL_VERSIONS 3.5)
    if( PYTHON_VERSION VERSION_LESS 3 )
        message(FATAL_ERROR "NGSolve supports only Python 3")
    endif( PYTHON_VERSION VERSION_LESS 3 )
    find_package(PythonInterp ${PYTHON_VERSION} REQUIRED)
    find_package(PythonLibs ${PYTHON_VERSION}  REQUIRED)

    set(PYTHON_LIBS "${PYTHON_LIBRARIES}")
    execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "from distutils.sysconfig import get_python_lib; print(get_python_lib(1,0,''))" OUTPUT_VARIABLE PYTHON_PACKAGES_INSTALL_DIR OUTPUT_STRIP_TRAILING_WHITESPACE)
    set_vars(PYTHON_LIBS PYTHON_PACKAGES_INSTALL_DIR PYBIND_INCLUDE_DIR PYTHON_INCLUDE_DIRS PYTHON_LIBRARIES PYTHON_EXECUTABLE PYTHON_VERSION)
endif (USE_PYTHON)

#######################################################################

if(USE_UMFPACK)
  ExternalProject_Add(
    suitesparse
    DEPENDS ${LAPACK_DEPENDENCIES}
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/umfpack
    GIT_REPOSITORY https://github.com/jlblancoc/suitesparse-metis-for-windows.git
    CMAKE_ARGS -DSUITESPARSE_USE_CUSTOM_BLAS_LAPACK_LIBS=ON -DSHARED=ON -DBUILD_METIS=OFF -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR} -DSUITESPARSE_CUSTOM_LAPACK_LIB=${LAPACK_LIBRARIES} -DSUITESPARSE_CUSTOM_BLAS_LIB=${LAPACK_LIBRARIES}
    )
  list(APPEND DEPENDENCIES suitesparse)
endif(USE_UMFPACK)

#######################################################################

if(USE_OCC AND WIN32)
    ExternalProject_Add(win_download_occ
      PREFIX ${CMAKE_CURRENT_BINARY_DIR}/tcl
      URL "http://www.asc.tuwien.ac.at/~mhochsteger/ngsuite/occ64.zip"
      UPDATE_COMMAND "" # Disable update
      BUILD_IN_SOURCE 1
      CONFIGURE_COMMAND ""
      BUILD_COMMAND ""
      INSTALL_COMMAND ${CMAKE_COMMAND} -E copy_directory . ${INSTALL_DIR}
      )
  list(APPEND DEPENDENCIES win_download_occ)
endif(USE_OCC AND WIN32)

#######################################################################

if(USE_GUI)
  include(cmake_modules/ExternalProject_TCLTK.cmake)
endif(USE_GUI)

#######################################################################
# propagate cmake variables to NGSolve subproject
set_vars(
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

ExternalProject_Add (ngsolve
  DEPENDS ${DEPENDENCIES}
  SOURCE_DIR ${PROJECT_SOURCE_DIR}
  CMAKE_ARGS -DUSE_SUPERBUILD=OFF ${EXTRA_CMAKE_ARGS} -DCMAKE_PREFIX_PATH=${INSTALL_DIR}
  INSTALL_COMMAND ""
  BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/ngsolve
  STEP_TARGETS build
)

install(CODE "execute_process(COMMAND cmake --build . --target install WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/ngsolve)")

