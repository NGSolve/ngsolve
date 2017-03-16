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
if(NOT WIN32 AND NOT INTEL_MIC)
    option( USE_NATIVE_ARCH  "build which -march=native" ON)
endif(NOT WIN32 AND NOT INTEL_MIC)

set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH}" "${CMAKE_CURRENT_SOURCE_DIR}/cmake/cmake_modules")
set(NETGEN_DIR "" CACHE PATH "Path to Netgen, leave empty to build Netgen automatically")

find_program(GIT_FOUND git)
if(GIT_FOUND AND EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/.git)
    execute_process(COMMAND git describe --tags --long --dirty WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} OUTPUT_VARIABLE git_version_string)
else()
  get_filename_component(git_version_string ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  string(REGEX REPLACE "^ngsolve_(.*)" "\\1" git_version_string "${git_version_string}")
endif()

string(REGEX REPLACE "^v([0-9]+)\\..*" "\\1" NGSOLVE_VERSION_MAJOR "${git_version_string}")
string(REGEX REPLACE "^v[0-9]+\\.([0-9]+).*" "\\1" NGSOLVE_VERSION_MINOR "${git_version_string}")
string(REGEX REPLACE "^v[0-9]+\\.[0-9]+\\.([0-9]+).*" "\\1" NGSOLVE_VERSION_PATCH "${git_version_string}")
string(REGEX REPLACE "^v[0-9]+\\.[0-9]+\\.[0-9]+\\-([0-9]+).*" "\\1" NGSOLVE_VERSION_TWEAK "${git_version_string}")
string(REGEX REPLACE "^v[0-9]+\\.[0-9]+\\.[0-9]+\\-[0-9]+\\-([0-9a-z]+).*" "\\1" NGSOLVE_VERSION_HASH "${git_version_string}")

set(NGSOLVE_VERSION_SHORT ${NGSOLVE_VERSION_MAJOR}.${NGSOLVE_VERSION_MINOR}.${NGSOLVE_VERSION_PATCH})
set(NGSOLVE_VERSION_LONG ${NGSOLVE_VERSION_SHORT}-${NGSOLVE_VERSION_TWEAK}-${NGSOLVE_VERSION_HASH})

if(NGSOLVE_VERSION_TWEAK)
  # no release version - nightly build
  set(NGSOLVE_VERSION ${NGSOLVE_VERSION_LONG})
else()
  # TWEAK is 0 -> current version has a tag assigned
  set(NGSOLVE_VERSION ${NGSOLVE_VERSION_SHORT})
endif()

set(NGSOLVE_VERSION_LONG ${NGSOLVE_VERSION_SHORT}-${NGSOLVE_VERSION_TWEAK}-${NGSOLVE_VERSION_HASH})

if (USE_SUPERBUILD)
  project (NGSuite)
  # execute the superbuild (this script will be invoked again without the
  # USE_SUPERBUILD option this time)
  include (cmake/SuperBuild.cmake)
  include (cmake/package.cmake)
  include(CPack)
  return()                      # stop processing this file further
else()
  project(NGSolve)
endif()

find_package(Netgen REQUIRED CONFIG HINTS ${NETGEN_DIR} ${NETGEN_DIR}/share/cmake ${INSTALL_DIR} ${INSTALL_DIR}/share/cmake ${INSTALL_DIR}/Contents/Resources/share/cmake)

# set(CPACK_DEBIAN_PACKAGE_VERSION "${NGSOLVE_VERSION_MAJOR}.${NGSOLVE_VERSION_MINOR}.${NGSOLVE_VERSION_PATCH}")

#######################################################################
include (CheckIncludeFiles)
check_include_files (dlfcn.h HAVE_DLFCN_H)
include (CheckCXXSymbolExists)
check_cxx_symbol_exists ( abi::__cxa_demangle "cxxabi.h" HAVE_CXA_DEMANGLE )
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/config.h.cmake ${CMAKE_CURRENT_BINARY_DIR}/config.h)

install (FILES ${CMAKE_CURRENT_BINARY_DIR}/config.h DESTINATION include COMPONENT ngsolve_devel )


#######################################################################

add_definitions(-DHAVE_NETGEN_SOURCES)
add_definitions(-DNGSOLVE_VERSION="${NGSOLVE_VERSION}")
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
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
endif(USE_NATIVE_ARCH)

#######################################################################
# append install paths of software in non-standard paths (e.g. openmpi, metis, intel mkl, ...)
# cmake -DUSE_MPI=ON -DCMAKE_PREFIX_PATH="/opt/openmpi165;/opt/metis51" ../
set(ADDITIONAL_PATHS "" CACHE PATH "List of paths to additional libraries in non-standard locations, separated by ';'")
if (ADDITIONAL_PATHS)
  set(CMAKE_PREFIX_PATH ${ADDITIONAL_PATHS})
endif (ADDITIONAL_PATHS)

#######################################################################
# use rpath
set(CMAKE_INSTALL_PREFIX ${NETGEN_DIR} CACHE INTERNAL "Prefix prepended to install directories" FORCE)
set(CMAKE_SKIP_BUILD_RPATH  FALSE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE) 
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
if("${isSystemDir}" STREQUAL "-1")
   set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
endif("${isSystemDir}" STREQUAL "-1")

#######################################################################
# build options
if(APPLE)
  set(ngs_install_dir_bin ../MacOS)
  set(ngs_install_dir_lib ../MacOS)
else(APPLE)
  set(ngs_install_dir_bin bin)
  set(ngs_install_dir_lib lib)
endif(APPLE)
set(ngs_install_dir EXPORT ngsolve-targets RUNTIME DESTINATION ${ngs_install_dir_bin} COMPONENT ngsolve LIBRARY DESTINATION ${ngs_install_dir_lib} COMPONENT ngsolve_devel ARCHIVE DESTINATION ${ngs_install_dir_lib} COMPONENT ngsolve_devel )
install(EXPORT ngsolve-targets DESTINATION share/cmake)

add_definitions(-DHAVE_CONFIG_H)
add_definitions(-DUSE_TIMEOFDAY)

include_directories (BEFORE ${PROJECT_SOURCE_DIR}/include)
include_directories (BEFORE ${PROJECT_BINARY_DIR})
include_directories (AFTER ${NETGEN_INCLUDE_DIRS})

link_directories ("${CMAKE_INSTALL_PREFIX}/${ngs_install_dir_lib}")

if(NETGEN_USE_GUI)
  install (FILES ngsolve.tcl DESTINATION ${ngs_install_dir_bin} COMPONENT ngsolve )
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
    add_definitions(-D_WIN32_WINNT=${ver} -DNOMINMAX)
    add_definitions(-DMSVC_EXPRESS -D_CRT_SECURE_NO_WARNINGS -DNGS_EXPORTS)

    # build convenience (aka object) libraries in windows)
    set(NGS_LIB_TYPE OBJECT)

    if(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
        # Intel Compiler
        set(CMAKE_CXX_FLAGS "/Od ${CMAKE_CXX_FLAGS} /MP /Qstd=c++14 /Qopt-report-file:report.log")
    else()
        # MS Visual Compiler
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP /bigobj")
    endif()

else(WIN32)
    # build shared libraries in windows)
    set(NGS_LIB_TYPE SHARED)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
endif(WIN32)
if(APPLE)
#    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -undefined dynamic_lookup")
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -undefined dynamic_lookup")
    add_definitions(-DMSG_NOSIGNAL=0)
#     set(CMAKE_SHARED_LIBRARY_SUFFIX ".so")
endif(APPLE)

#######################################################################
if (NETGEN_USE_GUI)
    add_definitions(-DTCL)
    include_directories(${NETGEN_TCL_INCLUDE_PATH})
endif (NETGEN_USE_GUI)

#######################################################################
if (NETGEN_USE_MPI)
    if(USE_MUMPS)
      enable_language(Fortran)
    endif(USE_MUMPS)
    find_package(MPI REQUIRED)
    add_definitions(-DPARALLEL)
    include_directories(${MPI_CXX_INCLUDE_PATH})
    set(CMAKE_EXE_LINKER_FLAGS    "${CMAKE_EXE_LINKER_FLAGS}    ${MPI_CXX_LINK_FLAGS} -ldl")
    set(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} ${MPI_CXX_LINK_FLAGS} -ldl")
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${MPI_CXX_LINK_FLAGS} -ldl")
endif (NETGEN_USE_MPI)
# set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wl,--no-undefined")

#######################################################################
if (USE_LAPACK)
    if(NOT LAPACK_LIBRARIES)
      find_package(LAPACK REQUIRED)
    endif()
    add_definitions(-DLAPACK)
endif (USE_LAPACK)

#######################################################################
if (USE_VT)
    add_definitions(-DVTRACE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -t:inst manual")
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
    add_definitions(-DUSE_MKL -DUSE_PARDISO -DLAPACK)


    if("${MKL_INTERFACE_LAYER}" STREQUAL "_ilp64")
        add_definitions(-DMKL_ILP64)
    endif()
endif (USE_MKL)

#######################################################################
if (USE_HYPRE)
  find_package(HYPRE REQUIRED)
  add_definitions(-DHYPRE)
  include_directories(${HYPRE_INCLUDES})
endif(USE_HYPRE)

if (USE_MUMPS)
    find_package(MUMPS REQUIRED)
    add_definitions(-DUSE_MUMPS)
    include_directories(${MUMPS_INCLUDE_DIR})
    link_directories(${MUMPS_DIR}/lib)
endif (USE_MUMPS)

#######################################################################
if (NETGEN_USE_PYTHON)
    add_definitions(-DNGS_PYTHON)
    include_directories(${NETGEN_PYBIND_INCLUDE_DIR})
    include_directories(${NETGEN_PYTHON_INCLUDE_DIRS})
endif (NETGEN_USE_PYTHON)

#######################################################################

if(USE_PARDISO)
  add_definitions(-DUSE_PARDISO)
  add_definitions(-DUSE_PARDISO400)
  if(NOT PARDISO_LIB)
      find_library(PARDISO_LIB NAMES libpardiso500-GNU481-X86-64.so)
  endif(NOT PARDISO_LIB)
  if(NOT PARDISO_LIB)
      message(FATAL_ERROR "Could NOT find Pardiso")
  endif(NOT PARDISO_LIB)
  find_package(OpenMP)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif(USE_PARDISO)

#######################################################################

if(USE_UMFPACK)
  find_package(SuiteSparse NO_MODULE HINTS ${UMFPACK_DIR})
  include(${USE_SuiteSparse})
  include_directories(${SuiteSparse_INCLUDE_DIRS})
  set(UMFPACK_LIBRARIES ${SuiteSparse_LIBRARIES} CACHE STRING "List of UMFPACK libraries")
  add_definitions(-DUSE_UMFPACK)
endif(USE_UMFPACK)

#######################################################################
if(USE_NUMA)
    add_definitions(-DUSE_NUMA)
    find_library(NUMA_LIB libnuma.so)
endif(USE_NUMA)


#######################################################################
if(USE_VTUNE)
    add_definitions(-DVTUNE)
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
get_directory_property( ngscxx_defines COMPILE_DEFINITIONS )
string(TOUPPER "CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE}" name)

foreach( d ${ngscxx_defines} )
    if(WIN32)
      if(NOT ${d} STREQUAL "NGS_EXPORTS")
        set( ngscxx_define_flags "${ngscxx_define_flags} /D${d}" )
      endif()
    else()
      set( ngscxx_define_flags "${ngscxx_define_flags} -D${d}" )
    endif()
endforeach()

set(ngscxx_flags "${${name}} ${CMAKE_CXX_FLAGS} -fPIC ${ngscxx_define_flags}")

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
    install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/pictures/logo.icns DESTINATION . RENAME Netgen.icns)
    if(CMAKE_OSX_SYSROOT)
        set( ngscxx_includes "-I${CMAKE_OSX_SYSROOT}/usr/include ${ngscxx_includes}" )
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
  set(ngsld_flags "${ngsld_flags} /LIBPATH:${INSTALL_DIR}/lib nglib.lib libngsolve.lib" )
else()
    set(ngsld_flags "${ngsld_flags} -L${INSTALL_DIR}/lib -Wl,--rpath=${INSTALL_DIR}/lib" )
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
    install (PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/ngscxx.bat DESTINATION ${ngs_install_dir_bin} COMPONENT ngsolve_devel )
    install (PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/ngsld.bat DESTINATION ${ngs_install_dir_bin} COMPONENT ngsolve_devel )
else(WIN32)
    get_property(MYCRULE GLOBAL PROPERTY RULE_LAUNCH_COMPILE)
    get_property(MYLRULE GLOBAL PROPERTY RULE_LAUNCH_LINK)
    file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/ngscxx "${ngscxx_header}
${CMAKE_CXX_COMPILER} ${ngscxx_flags} ${ngscxx_includes} $*\n")
    file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/ngsld "${ngsld_header}
${CMAKE_CXX_COMPILER} ${ngsld_flags} $*\n")
    install (PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/ngscxx DESTINATION ${ngs_install_dir_bin} COMPONENT ngsolve_devel )
    install (PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/ngsld DESTINATION ${ngs_install_dir_bin} COMPONENT ngsolve_devel )
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
    install (PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/ngspy DESTINATION ${ngs_install_dir_bin} COMPONENT ngsolve_devel )
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
  ${CMAKE_CURRENT_BINARY_DIR}/NGSolveConfig.cmake @ONLY)

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/NGSolveConfig.cmake DESTINATION share/cmake)


#######################################################################
# Configure message
# TODO: other message in case of failure
string(TOUPPER "CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE}" name)
set(flags "${${name}} ${CMAKE_CXX_FLAGS}")
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


