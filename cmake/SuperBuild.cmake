include (ExternalProject)

set_property (DIRECTORY PROPERTY EP_BASE dependencies)

set (DEPENDENCIES)
set (LAPACK_PROJECTS)
set (NETGEN_CMAKE_ARGS "" CACHE INTERNAL "")
set (NGSOLVE_CMAKE_ARGS "" CACHE INTERNAL "")
set (SUBPROJECT_CMAKE_ARGS "" CACHE INTERNAL "")

set (SUBPROJECT_ARGS
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/dependencies
)

# only show output on failure in ci-builds
if(DEFINED ENV{CI})
    set (SUBPROJECT_ARGS
        LOG_DOWNLOAD ON
        LOG_BUILD ON
        LOG_INSTALL ON
        LOG_CONFIGURE ON
    )
    if(${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.14.0")
        set (SUBPROJECT_ARGS
            ${SUBPROJECT_ARGS}
            LOG_OUTPUT_ON_FAILURE ON
            LOG_MERGED_STDOUTERR ON
        )
    endif()
endif()

macro(set_vars VAR_OUT)
  foreach(varname ${ARGN})
    if(NOT "${${varname}}" STREQUAL "")
      string(REPLACE ";" "|" varvalue "${${varname}}" )
      set(${VAR_OUT} ${${VAR_OUT}};-D${varname}=${varvalue} CACHE INTERNAL "")
    endif()
  endforeach()
endmacro()

# propagate all variables set on the command line using cmake -DFOO=BAR
# to Netgen and NGSolve subprojects
get_cmake_property(CACHE_VARS CACHE_VARIABLES)
foreach(CACHE_VAR ${CACHE_VARS})
  get_property(CACHE_VAR_HELPSTRING CACHE ${CACHE_VAR} PROPERTY HELPSTRING)
  if(CACHE_VAR_HELPSTRING STREQUAL "No help, variable specified on the command line." AND NOT CACHE_VAR STREQUAL "CMAKE_OSX_ARCHITECTURES")
    get_property(CACHE_VAR_TYPE CACHE ${CACHE_VAR} PROPERTY TYPE)
    set_vars(NETGEN_CMAKE_ARGS ${CACHE_VAR})
    set_vars(NGSOLVE_CMAKE_ARGS ${CACHE_VAR})
  endif()
endforeach()

set(NETGEN_DIR CACHE PATH "Path where Netgen is already installed. Setting this variable will skip the Netgen build and override the setting of CMAKE_INSTALL_PREFIX")

macro(set_flags_vars OUTPUT_VAR )
  foreach(varname ${ARGN})
    set_vars(${OUTPUT_VAR} ${varname} ${varname}_RELEASE ${varname}_MINSIZEREL ${varname}_RELWITHDEBINFO ${varname}_DEBUG)
  endforeach()
endmacro()

if(${CMAKE_GENERATOR} STREQUAL "Unix Makefiles")
  set(COMMON_BUILD_COMMAND $(MAKE) --silent )
else()
  set(COMMON_BUILD_COMMAND ${CMAKE_COMMAND} --build . --config ${CMAKE_BUILD_TYPE} )
endif()
#######################################################################

set_vars(SUBPROJECT_CMAKE_ARGS
    CMAKE_OSX_DEPLOYMENT_TARGET
    CMAKE_OSX_SYSROOT
    CMAKE_C_COMPILER
    CMAKE_CXX_COMPILER
    CMAKE_BUILD_TYPE
)

set(SUBPROJECT_CMAKE_ARGS "${SUBPROJECT_CMAKE_ARGS};-DCMAKE_POSITION_INDEPENDENT_CODE=ON" CACHE INTERNAL "")

if(USE_CCACHE)
  find_program(CCACHE_FOUND NAMES ccache ccache.bat)
  if(CCACHE_FOUND)
      set(SUBPROJECT_CMAKE_ARGS "${SUBPROJECT_CMAKE_ARGS};-DCMAKE_CXX_COMPILER_LAUNCHER=${CCACHE_FOUND}" CACHE INTERNAL "")
  endif()
endif()

#######################################################################
if(NETGEN_DIR)
  message(STATUS "Looking for NetgenConfig.cmake...")
  find_package(Netgen REQUIRED CONFIG HINTS ${NETGEN_DIR} ${NETGEN_DIR}/lib/cmake ${NETGEN_DIR}/share/cmake $ENV{NETGENDIR}/../share/cmake)
  if(NOT BUILD_FOR_CONDA)
  set(CMAKE_INSTALL_PREFIX ${NETGEN_DIR} CACHE PATH "Installation directory" FORCE)
  endif()
  add_custom_target(netgen_project)
else(NETGEN_DIR)
  message(STATUS "Build Netgen from git submodule")
  execute_process(COMMAND ${CMAKE_COMMAND} -P cmake/check_submodules.cmake WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} )
  add_custom_target(check_submodules_start ALL ${CMAKE_COMMAND} -P cmake/check_submodules.cmake WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} )
  add_custom_target(check_submodules_stop ALL ${CMAKE_COMMAND} -P cmake/check_submodules.cmake WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} DEPENDS ngsolve)

  if (USE_MPI AND NOT METIS_DIR)
    # hard-coded path to self-built metis
    message(STATUS "PARMETIS/METIS will be built by superbuild")
    set (METIS_DIR ${CMAKE_CURRENT_BINARY_DIR}/dependencies/parmetis)
    set (SUPERBUILD_METIS 1)
  endif (USE_MPI AND NOT METIS_DIR)

  if (USE_HYPRE AND NOT USE_MPI)
      message(FATAL_ERROR "Hypre needs MPI to be enabled (-DUSE_MPI=ON)!")
  endif (USE_HYPRE AND NOT USE_MPI)

  # propagate netgen-specific settings to Netgen subproject
  set_vars( NETGEN_CMAKE_ARGS
    CMAKE_INSTALL_PREFIX
    INSTALL_DEPENDENCIES
    INSTALL_PROFILES
    INTEL_MIC
    USE_CCACHE
    USE_NATIVE_ARCH
    ENABLE_UNIT_TESTS
    BUILD_STUB_FILES
    METIS_DIR
    )

  if (USE_MPI AND SUPERBUILD_METIS)
    # unset metis again
    unset(METIS_DIR)
  endif(USE_MPI AND SUPERBUILD_METIS)

  set_flags_vars(NETGEN_CMAKE_ARGS CMAKE_CXX_FLAGS CMAKE_SHARED_LINKER_FLAGS CMAKE_LINKER_FLAGS)
  ExternalProject_Add (netgen_project
    SOURCE_DIR ${PROJECT_SOURCE_DIR}/external_dependencies/netgen
    BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/netgen
    CONFIGURE_COMMAND ""
    LIST_SEPARATOR |
    ${SUBPROJECT_ARGS}
    BUILD_COMMAND ${COMMON_BUILD_COMMAND}
    CMAKE_ARGS
        ${NETGEN_CMAKE_ARGS}
        ${SUBPROJECT_CMAKE_ARGS}
    INSTALL_COMMAND ""
  )

  add_custom_target(install_netgen ALL
    ${CMAKE_COMMAND} --build ${CMAKE_CURRENT_BINARY_DIR}/netgen --target install --config ${CMAKE_BUILD_TYPE}
    DEPENDS netgen_project
  )

  list(APPEND DEPENDENCIES install_netgen)

  message("\n\nConfigure Netgen from submodule...")
  execute_process(COMMAND ${CMAKE_COMMAND} -G${CMAKE_GENERATOR} ${NETGEN_CMAKE_ARGS} ${SUBPROJECT_CMAKE_ARGS} ${PROJECT_SOURCE_DIR}/external_dependencies/netgen WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/netgen)
endif(NETGEN_DIR)

#######################################################################

set(LAPACK_LIBRARIES CACHE INTERNAL "Lapack libraries")
if(USE_MKL)
    set(MKL_MULTI_THREADED ON CACHE BOOL "Use threaded MKL libs")

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
    set_vars(NGSOLVE_CMAKE_ARGS MKL_ROOT MKL_STATIC MKL_SDL MKL_MULTI_THREADED MKL_INTERFACE_LAYER)
endif(USE_MKL)

if (USE_LAPACK)
    if(NOT LAPACK_LIBRARIES)
      if(WIN32)
        ExternalProject_Add(win_download_lapack
          URL "https://github.com/NGSolve/ngsolve_dependencies/releases/download/v1.0.0/lapack64.zip"
          URL_MD5 635432b6b41f23177b9116d4323c978c
          ${SUBPROJECT_ARGS}
          UPDATE_COMMAND "" # Disable update
          BUILD_IN_SOURCE 1
          CONFIGURE_COMMAND ""
          BUILD_COMMAND ""
          INSTALL_COMMAND ${CMAKE_COMMAND} -E copy_directory . ${CMAKE_INSTALL_PREFIX}
          )
        set(LAPACK_LIBRARIES ${CMAKE_INSTALL_PREFIX}/lib/BLAS.lib)
        list(APPEND LAPACK_PROJECTS win_download_lapack)
      else(WIN32)
        if(APPLE)
          set(LAPACK_LIBRARIES /System/Library/Frameworks/Accelerate.framework )
        else(APPLE)
          find_package(LAPACK REQUIRED)
        endif(APPLE)
      endif(WIN32)
    endif()
    set_vars(NGSOLVE_CMAKE_ARGS LAPACK_LIBRARIES)
endif (USE_LAPACK)

#######################################################################
if(USE_UMFPACK)
  if(BUILD_UMFPACK)
    set(UMFPACK_DIR ${CMAKE_CURRENT_BINARY_DIR}/umfpack/install CACHE PATH "Temporary directory to build UMFPACK")
    set(UMFPACK_STATIC ON)
    string(REPLACE ";" "|" LAPACK_LIBRARIES_LIST "${LAPACK_LIBRARIES}")
    ExternalProject_Add(
      suitesparse
      DEPENDS ${LAPACK_PROJECTS}
      GIT_REPOSITORY https://github.com/jlblancoc/suitesparse-metis-for-windows.git
      GIT_TAG 1618fd16ea34be287c0fbc32789ca280012a9280
      LIST_SEPARATOR |
      ${SUBPROJECT_ARGS}
      CMAKE_ARGS
          -DCMAKE_INSTALL_PREFIX=${UMFPACK_DIR}
          -DSHARED=OFF
          -DBUILD_METIS=OFF
          -DLAPACK_FOUND=FALSE # Use blas/lapack found by ngsolve
          -DSUITESPARSE_CUSTOM_BLAS_LIB=${LAPACK_LIBRARIES_LIST}
          -DSUITESPARSE_CUSTOM_LAPACK_LIB=${LAPACK_LIBRARIES_LIST}
          -DSUITESPARSE_INSTALL_PREFIX=${UMFPACK_DIR}
          -DSUITESPARSE_USE_CUSTOM_BLAS_LAPACK_LIBS=ON
         ${SUBPROJECT_CMAKE_ARGS}
      UPDATE_COMMAND ""
      )
    list(APPEND DEPENDENCIES suitesparse)
  endif()
  set_vars( NGSOLVE_CMAKE_ARGS UMFPACK_DIR UMFPACK_STATIC )
endif(USE_UMFPACK)

#######################################################################
if(USE_HYPRE AND NOT HYPRE_DIR)
  include(${CMAKE_CURRENT_LIST_DIR}/external_projects/hypre.cmake)
endif(USE_HYPRE AND NOT HYPRE_DIR)

#######################################################################
if(SUPERBUILD_METIS)
  message(STATUS "configuring PARMETIS/METIS")
  include(${CMAKE_CURRENT_LIST_DIR}/external_projects/parmetis.cmake)
  # we have to build parmetis (which also builds the metis library needed by netgen) before netgen
  add_dependencies(netgen_project project_parmetis)
endif(SUPERBUILD_METIS)

#######################################################################
if(USE_MUMPS AND NOT MUMPS_DIR)
  message(STATUS "configuring MUMPS")
  include(${CMAKE_CURRENT_LIST_DIR}/external_projects/mumps.cmake)
endif(USE_MUMPS AND NOT MUMPS_DIR)

#######################################################################
# propagate cmake variables to NGSolve subproject
set_vars( NGSOLVE_CMAKE_ARGS
  CMAKE_INSTALL_PREFIX

  USE_LAPACK
  USE_MPI
  USE_CUDA
  USE_MKL
  USE_HYPRE
  USE_MUMPS
  USE_PARDISO
  USE_UMFPACK
  USE_VTUNE
  USE_CCACHE
  USE_NATIVE_ARCH
  NETGEN_DIR
  Netgen_DIR
  INSTALL_DEPENDENCIES 
  INTEL_MIC
  ENABLE_UNIT_TESTS
  BUILD_STUB_FILES
  )

set_flags_vars(NGSOLVE_CMAKE_ARGS CMAKE_CXX_FLAGS CMAKE_SHARED_LINKER_FLAGS CMAKE_LINKER_FLAGS)

set(NGSOLVE_PREFIX_PATH ${CMAKE_INSTALL_PREFIX} ${CMAKE_PREFIX_PATH})

string(REPLACE ";" "|" NGSOLVE_PREFIX_PATH_ALT_SEP "${NGSOLVE_PREFIX_PATH}")
ExternalProject_Add (ngsolve
  DEPENDS ${DEPENDENCIES} ${LAPACK_PROJECTS}
  SOURCE_DIR ${PROJECT_SOURCE_DIR}
  LIST_SEPARATOR |
  CMAKE_ARGS
         ${NGSOLVE_CMAKE_ARGS}
         -DUSE_SUPERBUILD=OFF
         -DCMAKE_PREFIX_PATH=${NGSOLVE_PREFIX_PATH}
         ${SUBPROJECT_CMAKE_ARGS}
  INSTALL_COMMAND ""
  BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/ngsolve
  BUILD_COMMAND ${COMMON_BUILD_COMMAND}
  ${SUBPROJECT_ARGS}
  )


install(CODE "execute_process(COMMAND \"${CMAKE_COMMAND}\" --build . --config ${CMAKE_BUILD_TYPE} --target install WORKING_DIRECTORY \"${CMAKE_CURRENT_BINARY_DIR}/ngsolve\")")

add_custom_target(test_ngsolve
  ${CMAKE_COMMAND} --build ${CMAKE_CURRENT_BINARY_DIR}/ngsolve
                   --target test
                   --config ${CMAKE_BUILD_TYPE}
                   )

# Check if the git submodules (i.e. netgen) are up to date
# in case, something is wrong, emit a warning but continue
 ExternalProject_Add_Step(ngsolve check_submodules
   COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_LIST_DIR}/check_submodules.cmake
   DEPENDERS install # Steps on which this step depends
   )

# Due to 'ALWAYS 1', this step is always run which also forces a build of
# the ngsolve subproject
 ExternalProject_Add_Step(ngsolve check_submodules1
   COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_LIST_DIR}/check_submodules.cmake
   DEPENDEES configure # Steps on which this step depends
   DEPENDERS build     # Steps that depend on this step
   ALWAYS 1            # No stamp file, step always runs
   )



if(WIN32)
  file(TO_NATIVE_PATH ${CMAKE_INSTALL_PREFIX}/bin netgendir)
  file(TO_NATIVE_PATH ${CMAKE_INSTALL_PREFIX}/${NETGEN_PYTHON_PACKAGES_INSTALL_DIR} pythonpath)
  add_custom_target(set_netgendir
    setx NETGENDIR  ${netgendir}
  )
  add_custom_target(set_pythonpath
    setx PYTHONPATH "${pythonpath};$ENV{PYTHONPATH}"
  )
  add_custom_target(set_environment_variables)
  add_dependencies(set_environment_variables set_netgendir set_pythonpath)
endif(WIN32)
 
