if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "RelWithDebInfo" CACHE STRING INTERNAL)
endif(NOT CMAKE_BUILD_TYPE)

cmake_minimum_required(VERSION 2.8.12)

# Handle superbuild first
option( USE_SUPERBUILD   "Whether or not a superbuild should be invoked" ON)
option( USE_LAPACK       "build ngsolve with LAPACK" ON)
option( USE_VT           "enable vampirtrace")
# option( USE_CUDA         "enable CUDA GPU support")
option( USE_MKL          "enable MKL")
option( USE_HYPRE        "enable HYPRE support")
option( USE_MUMPS        "enable sparse direct solver MUMPS")
option( USE_PARDISO      "enable pardiso sparse direct solver")
option( USE_UMFPACK      "enable umfpack sparse direct solver")
option( INTEL_MIC        "cross compile for intel xeon phi")
option( USE_VTUNE        "include vtune pause/resume numproc")
option( USE_NUMA         "compile with NUMA-aware code")
option( USE_CCACHE       "use ccache")
option( INSTALL_DEPENDENCIES "install dependencies like netgen or solver libs, useful for packaging" OFF )
option( ENABLE_UNIT_TESTS "Enable Catch unit tests")
if(NOT WIN32 AND NOT INTEL_MIC)
    option( USE_NATIVE_ARCH  "build which -march=native" ON)
endif(NOT WIN32 AND NOT INTEL_MIC)

set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH}" "${CMAKE_CURRENT_SOURCE_DIR}/cmake/cmake_modules")
set(NETGEN_DIR "" CACHE PATH "Path to Netgen, leave empty to build Netgen automatically")

include (${CMAKE_CURRENT_LIST_DIR}/cmake/generate_version_file.cmake)

if(APPLE)
  set(INSTALL_DIR_DEFAULT /Applications/Netgen.app)
else(APPLE)
  if(WIN32)
    set(INSTALL_DIR_DEFAULT "C:/netgen")
  else(WIN32)
    set(INSTALL_DIR_DEFAULT /opt/netgen)
  endif(WIN32)
endif(APPLE)

if(INSTALL_DIR)
  message(WARNING "INSTALL_DIR is deprecated, use CMAKE_INSTALL_PREFIX instead")
  set(INSTALL_DIR_DEFAULT ${INSTALL_DIR})
endif(INSTALL_DIR)

if (USE_SUPERBUILD)
  project (NGSuite)
  if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX "${INSTALL_DIR_DEFAULT}" CACHE PATH "Install directory" FORCE)
  endif()
  # execute the superbuild (this script will be invoked again without the
  # USE_SUPERBUILD option this time)
  include (cmake/SuperBuild.cmake)
  include (cmake/package.cmake)
  include(CPack)
  return()                      # stop processing this file further
else()
  project(NGSolve)
  if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX "${INSTALL_DIR_DEFAULT}" CACHE PATH "Install directory" FORCE)
  endif()
endif()

find_package(Netgen REQUIRED CONFIG HINTS ${CMAKE_INSTALL_PREFIX}/Contents/Resources)

set(NGSOLVE_INSTALL_DIR_PYTHON  ${NETGEN_INSTALL_DIR_PYTHON})
set(NGSOLVE_INSTALL_DIR_BIN     ${NETGEN_INSTALL_DIR_BIN})
set(NGSOLVE_INSTALL_DIR_LIB     ${NETGEN_INSTALL_DIR_LIB})
set(NGSOLVE_INSTALL_DIR_INCLUDE ${NETGEN_INSTALL_DIR_INCLUDE})
set(NGSOLVE_INSTALL_DIR_RES     ${NETGEN_INSTALL_DIR_RES})
string(REPLACE "netgen" "ngsolve" NGSOLVE_INSTALL_DIR_CMAKE ${NETGEN_INSTALL_DIR_CMAKE})

#######################################################################

set(NGSOLVE_COMPILE_INCLUDE_DIRS "")
set(NGSOLVE_COMPILE_DEFINITIONS "")
set(NGSOLVE_COMPILE_OPTIONS "")
list(APPEND NGSOLVE_COMPILE_DEFINITIONS HAVE_NETGEN_SOURCES)

#######################################################################

include (CheckIncludeFiles)
check_include_files (dlfcn.h HAVE_DLFCN_H)
include (CheckCXXSymbolExists)
check_cxx_symbol_exists ( abi::__cxa_demangle "cxxabi.h" HAVE_CXA_DEMANGLE )
if(HAVE_DLFCN_H)
    list(APPEND NGSOLVE_COMPILE_DEFINITIONS HAVE_DLFCN_H)
endif()
if(HAVE_CXA_DEMANGLE)
    list(APPEND NGSOLVE_COMPILE_DEFINITIONS HAVE_CXA_DEMANGLE)
endif()

add_custom_target(generate_version_file
  cmake
    -DBDIR=${CMAKE_CURRENT_BINARY_DIR}
    -P ${CMAKE_CURRENT_LIST_DIR}/cmake/generate_version_file.cmake
    )

#######################################################################
if(USE_CCACHE)
  find_program(CCACHE_FOUND NAMES ccache ccache.bat)
  if(CCACHE_FOUND)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE ${CCACHE_FOUND})
  endif(CCACHE_FOUND)
endif(USE_CCACHE)

#######################################################################
if(INTEL_MIC)
    set(MKL_ARCH "mic")
    include(cmake/mic.cmake)
endif(INTEL_MIC)

#######################################################################
if(USE_NATIVE_ARCH)
    list(APPEND NGSOLVE_COMPILE_OPTIONS -march=native)
endif(USE_NATIVE_ARCH)

#######################################################################
if(ENABLE_UNIT_TESTS)
  include(${CMAKE_CURRENT_LIST_DIR}/cmake/external_projects/catch.cmake)
endif(ENABLE_UNIT_TESTS)

#######################################################################
# append install paths of software in non-standard paths (e.g. openmpi, metis, intel mkl, ...)
# cmake -DUSE_MPI=ON -DCMAKE_PREFIX_PATH="/opt/openmpi165;/opt/metis51" ../
set(ADDITIONAL_PATHS "" CACHE PATH "List of paths to additional libraries in non-standard locations, separated by ';'")
if (ADDITIONAL_PATHS)
  set(CMAKE_PREFIX_PATH ${ADDITIONAL_PATHS})
endif (ADDITIONAL_PATHS)

#######################################################################
# use rpath
set(CMAKE_MACOSX_RPATH TRUE)
set(CMAKE_INSTALL_RPATH "${NETGEN_RPATH_TOKEN};${NETGEN_RPATH_TOKEN}/${NETGEN_RPATH}")

#######################################################################
set(ngs_install_dir EXPORT ngsolve-targets RUNTIME DESTINATION ${NGSOLVE_INSTALL_DIR_BIN} COMPONENT ngsolve LIBRARY DESTINATION ${NGSOLVE_INSTALL_DIR_LIB} COMPONENT ngsolve_devel ARCHIVE DESTINATION ${NGSOLVE_INSTALL_DIR_LIB} COMPONENT ngsolve_devel )
install(EXPORT ngsolve-targets DESTINATION ${NGSOLVE_INSTALL_DIR_CMAKE})

list(APPEND NGSOLVE_COMPILE_DEFINITIONS USE_TIMEOFDAY)
list(APPEND NGSOLVE_INCLUDE_DIRS $<INSTALL_INTERFACE:${NETGEN_INSTALL_DIR_INCLUDE}>)
list(APPEND NGSOLVE_INCLUDE_DIRS $<INSTALL_INTERFACE:${NETGEN_INSTALL_DIR_INCLUDE}/include>)

include_directories (BEFORE ${PROJECT_SOURCE_DIR}/include)
include_directories (BEFORE ${PROJECT_BINARY_DIR})
include_directories (AFTER ${NETGEN_INCLUDE_DIRS})

if(NETGEN_USE_GUI)
  file(RELATIVE_PATH BIN_TO_LIB_RELPATH ${NETGEN_BINARY_DIR} ${NETGEN_LIBRARY_DIR})
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/ngsolve.tcl
    ${CMAKE_CURRENT_BINARY_DIR}/ngsolve.tcl @ONLY)

  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/ngsolve.tcl DESTINATION ${NGSOLVE_INSTALL_DIR_BIN})
endif(NETGEN_USE_GUI)

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DDEBUG -Wall")

#######################################################################
# platform specific definitions and flags
macro(get_WIN32_WINNT version)
    if (WIN32 AND CMAKE_SYSTEM_VERSION)
        if("${CMAKE_SYSTEM_VERSION}" MATCHES "^([0-9]+)\\.([0-9]+)")
            math(EXPR ver "${CMAKE_MATCH_1}*100 + ${CMAKE_MATCH_2}")
        endif()
        set(${version} "0x${ver}")
    endif()
endmacro()
get_WIN32_WINNT(ver)
if(WIN32)
    add_definitions(-DNGS_EXPORTS)
    list(APPEND NGSOLVE_COMPILE_DEFINITIONS NOMINMAX _WIN32_WINNT=${ver} MSVC_EXPRESS _CRT_SECURE_NO_WARNINGS )

    # build convenience (aka object) libraries in windows)
    set(NGS_LIB_TYPE OBJECT)

    if(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
        # Intel Compiler
        list(APPEND NGSOLVE_COMPILE_OPTIONS /Od /MP /Qstd=c++14 /Qopt-report-file:report.log)
    else()
        # MS Visual Compiler
        list(APPEND NGSOLVE_COMPILE_OPTIONS /MP /bigobj)
    endif()

else(WIN32)
    # build shared libraries in windows)
    set(NGS_LIB_TYPE SHARED)
    list(APPEND NGSOLVE_COMPILE_OPTIONS -std=c++14)
endif(WIN32)
if(APPLE)
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -undefined dynamic_lookup")
    list(APPEND NGSOLVE_COMPILE_DEFINITIONS MSG_NOSIGNAL=0)
endif(APPLE)

#######################################################################
if (NETGEN_USE_GUI)
    list(APPEND NGSOLVE_COMPILE_DEFINITIONS TCL)
    list(APPEND NGSOLVE_INCLUDE_DIRS ${NETGEN_TCL_INCLUDE_PATH})
endif (NETGEN_USE_GUI)

#######################################################################
if (NETGEN_USE_MPI)
    if(USE_MUMPS)
      enable_language(Fortran)
    endif(USE_MUMPS)
    find_package(MPI REQUIRED)
    list(APPEND NGSOLVE_COMPILE_DEFINITIONS PARALLEL)
    include_directories(${MPI_CXX_INCLUDE_PATH})
    set(CMAKE_EXE_LINKER_FLAGS    "${CMAKE_EXE_LINKER_FLAGS}    ${MPI_CXX_LINK_FLAGS} -ldl")
    set(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} ${MPI_CXX_LINK_FLAGS} -ldl")
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${MPI_CXX_LINK_FLAGS} -ldl")
endif (NETGEN_USE_MPI)

#######################################################################
if (USE_LAPACK)
    if(NOT LAPACK_LIBRARIES)
      find_package(LAPACK REQUIRED)
    endif()
    list(APPEND NGSOLVE_COMPILE_DEFINITIONS LAPACK)
endif (USE_LAPACK)

#######################################################################
if (USE_VT)
    list(APPEND NGSOLVE_COMPILE_DEFINITIONS VTRACE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ")
    list(APPEND NGSOLVE_COMPILE_OPTIONS -t:inst manual)
endif (USE_VT)

#######################################################################
if (USE_CUDA)
    find_package(CUDA)
    set(CUDA_SEPARABLE_COMPILATION ON)
    set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS};-std=c++14;-O2;-DUSE_MYCOMPLEX;-DHD="__host__ __device__")
    set(CUDA_OPTIONS -Xcompiler -fPIC)
    SET (EXTRA_LIBS ${EXTRA_LIBS} cudart cublas cusparse)
endif (USE_CUDA)

#######################################################################
if (USE_MKL)
    find_package(MKL REQUIRED)
    include_directories(${MKL_INCLUDE_DIRS})
    list(APPEND NGSOLVE_COMPILE_DEFINITIONS USE_MKL USE_PARDISO LAPACK)


    if("${MKL_INTERFACE_LAYER}" STREQUAL "_ilp64")
      list(APPEND NGSOLVE_COMPILE_DEFINITIONS MKL_ILP64)
    endif()
endif (USE_MKL)

#######################################################################
if (USE_HYPRE)
  find_package(HYPRE REQUIRED)
  list(APPEND NGSOLVE_COMPILE_DEFINITIONS HYPRE)
  include_directories(${HYPRE_INCLUDES})
endif(USE_HYPRE)

if (USE_MUMPS)
    find_package(MUMPS REQUIRED)
    list(APPEND NGSOLVE_COMPILE_DEFINITIONS USE_MUMPS)
    include_directories(${MUMPS_INCLUDE_DIR})
    link_directories(${MUMPS_DIR}/lib)
endif (USE_MUMPS)

#######################################################################
if (NETGEN_USE_PYTHON)
    list(APPEND NGSOLVE_COMPILE_DEFINITIONS NGS_PYTHON)
    list(APPEND NGSOLVE_INCLUDE_DIRS ${NETGEN_PYBIND_INCLUDE_DIR} ${NETGEN_PYTHON_INCLUDE_DIRS})
endif (NETGEN_USE_PYTHON)

#######################################################################

if(USE_PARDISO)
  list(APPEND NGSOLVE_COMPILE_DEFINITIONS USE_PARDISO)
  list(APPEND NGSOLVE_COMPILE_DEFINITIONS USE_PARDISO400)
  if(NOT PARDISO_LIB)
      find_library(PARDISO_LIB NAMES libpardiso500-GNU481-X86-64.so)
  endif(NOT PARDISO_LIB)
  if(NOT PARDISO_LIB)
      message(FATAL_ERROR "Could NOT find Pardiso")
  endif(NOT PARDISO_LIB)
  find_package(OpenMP)
  list(APPEND NGSOLVE_COMPILE_OPTIONS ${OpenMP_CXX_FLAGS})
endif(USE_PARDISO)

#######################################################################

if(USE_UMFPACK)
  find_package(SuiteSparse NO_MODULE HINTS ${UMFPACK_DIR})
  include(${USE_SuiteSparse})
  include_directories(${SuiteSparse_INCLUDE_DIRS})
  set(UMFPACK_LIBRARIES ${SuiteSparse_LIBRARIES} CACHE STRING "List of UMFPACK libraries")
  list(APPEND NGSOLVE_COMPILE_DEFINITIONS USE_UMFPACK)
endif(USE_UMFPACK)

#######################################################################
if(USE_NUMA)
    list(APPEND NGSOLVE_COMPILE_DEFINITIONS USE_NUMA)
    find_library(NUMA_LIB libnuma.so)
endif(USE_NUMA)


#######################################################################
if(USE_VTUNE)
    list(APPEND NGSOLVE_COMPILE_DEFINITIONS VTUNE)
    find_path(VTUNE_INCLUDE_DIR libittnotify.h)
    include_directories(${VTUNE_INCLUDE_DIR})
    get_filename_component(VTUNE_DIR ${VTUNE_INCLUDE_DIR} PATH)
    if(INTEL_MIC)
        find_library(VTUNE_LIB libittnotify.a PATH ${VTUNE_DIR}/bin64/k1om)
    else(INTEL_MIC)
        find_library(VTUNE_LIB libittnotify.a PATH ${VTUNE_DIR}/lib64)
    endif(INTEL_MIC)
endif(USE_VTUNE)

#######################################################################
if(${CMAKE_BUILD_TYPE} STREQUAL "Debug" OR ${CMAKE_BUILD_TYPE} STREQUAL "DEBUG")
  set(NGS_TEST_TIMEOUT 300)
else()
  set(NGS_TEST_TIMEOUT 60)
endif()

#######################################################################

add_subdirectory(include)
add_subdirectory(ngstd)
add_subdirectory(basiclinalg)
add_subdirectory(parallel)
add_subdirectory(linalg)
add_subdirectory(fem)
add_subdirectory(multigrid)
add_subdirectory(comp)
add_subdirectory(solve)
add_subdirectory(python)
add_subdirectory(pde_tutorial)
add_subdirectory(py_tutorials)
add_subdirectory(tests)
add_subdirectory(docs)

#######################################################################
# generate ngscxx compiler wrapper

# compile flags and preprocessor definitions
string(TOUPPER "CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE}" name)

foreach( d ${NGSOLVE_COMPILE_DEFINITIONS} )
    if(WIN32)
      if(NOT ${d} STREQUAL "NGS_EXPORTS")
        set( ngscxx_define_flags "${ngscxx_define_flags} /D${d}" )
      endif()
    else()
      set( ngscxx_define_flags "${ngscxx_define_flags} -D${d}" )
    endif()
endforeach()

set(ngscxx_flags "${${name}} ${CMAKE_CXX_FLAGS} ${ngscxx_define_flags}")

foreach( opt ${NGSOLVE_COMPILE_OPTIONS} )
    set(ngscxx_flags "${ngscxx_flags} ${opt}")
endforeach()
if(NOT WIN32)
  set(ngscxx_flags "${ngscxx_flags} -fPIC")
endif()

# include directories
set(ngscxx_include_dirs ${NETGEN_INCLUDE_DIRS} ${CMAKE_INSTALL_PREFIX}/include ${NETGEN_TCL_INCLUDE_PATH} ${NETGEN_MPI_CXX_INCLUDE_PATH} ${MKL_INCLUDE_DIRS} ${MUMPS_INCLUDE_DIR} ${NETGEN_PYTHON_INCLUDE_DIRS} ${UMFPACK_INCLUDE_DIR} ${VTUNE_INCLUDE_DIR})


foreach( d ${ngscxx_include_dirs} )
  if(WIN32)
    set( ngscxx_includes "${ngscxx_includes} /I\"${d}\"" )
  else()
    set( ngscxx_includes "${ngscxx_includes} -I${d}" )
  endif()
endforeach()

# linker flags
if(USE_MKL AND NOT INSTALL_DEPENDENCIES)
    set( ngsld_flags "${ngsld_flags} ${MKL_LINK_FLAGS}" )
endif()

if(APPLE)
    install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/pictures/logo.icns DESTINATION ${NGSOLVE_INSTALL_DIR_RES}/../ RENAME Netgen.icns)
    if(CMAKE_OSX_SYSROOT)
        set( ngscxx_includes "-isysroot ${CMAKE_OSX_SYSROOT} ${ngscxx_includes}" )
    endif()
    set(ngscxx_flags "${ngscxx_flags} -Qunused-arguments")
    set(ngscxx_includes "${ngscxx_includes} -I$Netgen_BUNDLE/Contents/Resources/include" )  
    set(ngsld_flags "${ngsld_flags} -L$Netgen_BUNDLE/Contents/MacOS  -undefined dynamic_lookup" )
    set(ngscxx_header "\
#!/bin/sh
Netgen_BUNDLE=\"`echo \"$0\" | sed -e 's/\\/Contents\\/MacOS\\/ngscxx//'`\"
Netgen_MACOS=\"$Netgen_BUNDLE/Contents/MacOS\"
")
    set(ngsld_header "\
#!/bin/sh
Netgen_BUNDLE=\"`echo \"$0\" | sed -e 's/\\/Contents\\/MacOS\\/ngsld//'`\"
Netgen_MACOS=\"$Netgen_BUNDLE/Contents/MacOS\"
")
elseif(WIN32)
  set(ngsld_flags "${ngsld_flags} /LIBPATH:${NETGEN_LIBRARY_DIR} nglib.lib libngsolve.lib" )
else()
    set(ngsld_flags "${ngsld_flags} -L${NETGEN_LIBRARY_DIR} -Wl,--rpath=${NETGEN_LIBRARY_DIR}" )
    set(ngscxx_header "#!/bin/sh")
    set(ngsld_header "#!/bin/sh")
endif()

if(WIN32)
    get_filename_component(VC_BIN_DIR ${CMAKE_CXX_COMPILER} PATH)
    get_filename_component(VC_DIR ${VC_BIN_DIR}/../../ ABSOLUTE)
    file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/ngscxx.bat 
      "call \"${VC_DIR}/vcvarsall.bat\" amd64\ncl /c ${ngscxx_flags} ${ngscxx_includes} %*\n")
    file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/ngsld.bat
      "call \"${VC_DIR}/vcvarsall.bat\" amd64\nlink /DLL %* ${ngsld_flags}\n")
    install (PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/ngscxx.bat DESTINATION ${NGSOLVE_INSTALL_DIR_BIN} COMPONENT ngsolve_devel )
    install (PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/ngsld.bat DESTINATION ${NGSOLVE_INSTALL_DIR_BIN} COMPONENT ngsolve_devel )
else(WIN32)
    set(ngscxx_compile_command "${CMAKE_CXX_COMPILER} ${ngscxx_flags} ${ngscxx_includes} $*")
    if(CCACHE_FOUND)
      set(ngscxx_compile_command "
if [ -f ${CCACHE_FOUND} ]
  then ${CCACHE_FOUND} ${ngscxx_compile_command}
  else ${ngscxx_compile_command}
fi")
    endif(CCACHE_FOUND)
    file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/ngscxx "${ngscxx_header}
${ngscxx_compile_command}\n")
    file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/ngsld "${ngsld_header}
${CMAKE_CXX_COMPILER} ${ngsld_flags} $*\n")
install (PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/ngscxx DESTINATION ${NGSOLVE_INSTALL_DIR_BIN} COMPONENT ngsolve_devel )
install (PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/ngsld DESTINATION ${NGSOLVE_INSTALL_DIR_BIN} COMPONENT ngsolve_devel )
endif(WIN32)

if (NETGEN_USE_MPI)
    # generate python3 wrapper for mpi (to preload mkl and mpi libs)
    set(ngspy_preload_libs "${MPI_C_LIBRARIES}")
    if(USE_MKL AND NOT MKL_STATIC)
      set(ngspy_preload_libs "${ngspy_preload_libs}:${MKL_CORE_LIBRARY}:${MKL_THREADING_LIBRARY}:${MKL_INTERFACE_LIBRARY}:${MKL_BLACS_LIBRARY}")
    endif()
    file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/ngspy
"#! /usr/bin/bash
LD_PRELOAD=$LD_PRELOAD:${ngspy_preload_libs} ${NETGEN_PYTHON_EXECUTABLE} $*
")
    install (PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/ngspy DESTINATION ${NGSOLVE_INSTALL_DIR_BIN} COMPONENT ngsolve_devel )
endif (NETGEN_USE_MPI)




#######################################################################
# CTest
enable_testing()
include(CTest)
add_custom_target(performance_tests ctest -L performance -D Experimental)


#######################################################################
# uninstall target
configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/cmake_uninstall.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
    IMMEDIATE @ONLY)
# add_custom_target(uninstall
#     COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cmake/NGSolveConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/NGSolveConfig.cmake @ONLY ESCAPE_QUOTES)

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/NGSolveConfig.cmake DESTINATION ${NGSOLVE_INSTALL_DIR_CMAKE})


#######################################################################
# Configure message
# TODO: other message in case of failure
string(TOUPPER "CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE}" name)
set(flags "${ngscxx_flags}")
message("
------------------------------------------------------------------------
  ${PROJECT_NAME} ${NGSOLVE_VERSION}:  Automatic configuration OK.

  Build type: ${CMAKE_BUILD_TYPE}
  Compiler: ${CMAKE_CXX_COMPILER}
  Flags: ${flags}

  Enabled functionality:

    GUI: ............... ${NETGEN_USE_GUI}
    LAPACK: ............ ${USE_LAPACK}
    MKL: ............... ${USE_MKL}
    PYTHON: ............ ${NETGEN_USE_PYTHON}
    MPI: ............... ${NETGEN_USE_MPI}
    CUDA: .............. ${USE_CUDA}
    MUMPS: ............. ${USE_MUMPS}
    PARDISO: ........... ${USE_PARDISO}
    INTEL_MIC: ......... ${INTEL_MIC}
    VTUNE: ............. ${USE_VTUNE}
    VAMPIRTRACE: ....... ${USE_VT}


  Building:
  ")
if(WIN32)
message("
    Open ${CMAKE_BINARY_DIR}/NGSolve.sln and build solution to compile ${PROJECT_NAME}.

    Build \"INSTALL\" to install ${PROJECT_NAME}.
")
else(WIN32)
message("
    Type 'make' to compile ${PROJECT_NAME}.

    Type 'make install' to install ${PROJECT_NAME}.

")
endif(WIN32)
message("
  Install directory:

    ${CMAKE_INSTALL_PREFIX}")
message("
------------------------------------------------------------------------
")

if (ADDITIONAL_PATHS)
  message(WARNING "The use of ADDITIONAL_PATHS is depricated, use CMAKE_PREFIX_PATH instead.")
endif (ADDITIONAL_PATHS)


