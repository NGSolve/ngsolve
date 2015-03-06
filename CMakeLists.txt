if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "RelWithDebInfo" CACHE STRING INTERNAL)
endif(NOT CMAKE_BUILD_TYPE)

cmake_minimum_required(VERSION 2.8)
project(NGSolve)
set(NGSOLVE_VERSION_MAJOR 6)
set(NGSOLVE_VERSION_MINOR 1)
set(NGSOLVE_VERSION "${NGSOLVE_VERSION_MAJOR}.${NGSOLVE_VERSION_MINOR}-dev")
set(PACKAGE_VERSION ${NGSOLVE_VERSION})

if(NOT WIN32)
    option( USE_NATIVE_ARCH  "build which -march=native" ON)
endif(NOT WIN32)

option( USE_GUI          "build ngsolve with GUI" ON)
option( USE_PYTHON       "enable python shell" ON)
option( USE_LAPACK       "build ngsolve with LAPACK" ON)
option( USE_MPI          "enable mpi parallelization")
option( USE_VT           "enable vampirtrace")
option( USE_CUDA         "enable CUDA GPU support")
option( USE_MKL          "enable MKL")
option( USE_MUMPS        "enable sparse direct solver MUMPS")
option( USE_GOLD         "enable gold")
option( INTEL_MIC        "cross compile for intel xeon phi")
option( USE_VTUNE        "include vtune pause/resume numproc")

#######################################################################
if(USE_NATIVE_ARCH)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
endif(USE_NATIVE_ARCH)

#######################################################################
# find netgen
if(NOT INSTALL_DIR)
    set(netgen_env $ENV{NETGENDIR})
    if(NOT netgen_env)
        message( FATAL_ERROR "Could not find Netgen. Set either NETGENDIR environment variable or INSTALL_DIR in the CMake cache.")
    else(NOT netgen_env)
        get_filename_component(NETGEN_BASEDIR $ENV{NETGENDIR} PATH)
        set(INSTALL_DIR ${NETGEN_BASEDIR} CACHE PATH "Path to install NGSolve (i.e. where netgen is installed)")
    endif(NOT netgen_env)
endif(NOT INSTALL_DIR)
set(CMAKE_INSTALL_PREFIX "${INSTALL_DIR}" CACHE INTERNAL "Prefix prepended to install directories" FORCE)

#######################################################################
# append install paths of software in non-standard paths (e.g. openmpi, metis, intel mkl, ...)
# cmake -DUSE_MPI=ON -DADDITIONAL_PATHS="/opt/openmpi165;/opt/metis51" ../
set(ADDITIONAL_PATHS "" CACHE PATH "List of paths to additional libraries in non-standard locations, seperated by ';'")
if (ADDITIONAL_PATHS)
    list(REVERSE CMAKE_SYSTEM_PREFIX_PATH)
    list(APPEND CMAKE_SYSTEM_PREFIX_PATH ${ADDITIONAL_PATHS})
    list(REVERSE CMAKE_SYSTEM_PREFIX_PATH)
endif (ADDITIONAL_PATHS)

#######################################################################
# use rpath
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
set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH}" "${PROJECT_SOURCE_DIR}/cmake_modules")
set(ngs_install_dir RUNTIME DESTINATION bin LIBRARY DESTINATION lib ARCHIVE DESTINATION lib)

add_definitions(-DHAVE_CONFIG_H)
add_definitions(-DUSE_TIMEOFDAY)

set(CMAKE_INCLUDE_CURRENT_DIR ON)
include_directories ("${PROJECT_SOURCE_DIR}/include")
include_directories ("${PROJECT_BINARY_DIR}")
include_directories ("${INSTALL_DIR}/include")

link_directories ("${INSTALL_DIR}/lib")

find_package(OpenMP)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DDEBUG")

include (CheckIncludeFiles)
check_include_files (dlfcn.h HAVE_DLFCN_H)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/config.h.cmake ${CMAKE_CURRENT_BINARY_DIR}/config.h)

install (FILES ${CMAKE_CURRENT_BINARY_DIR}/config.h DESTINATION include)
install (FILES ngsolve.tcl DESTINATION bin)
install (PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/ngscxx DESTINATION bin)

#######################################################################
# platform specific definitions and flags
macro(get_WIN32_WINNT version)
    if (WIN32 AND CMAKE_SYSTEM_VERSION)
        set(ver ${CMAKE_SYSTEM_VERSION})
        string(REPLACE "." "" ver ${ver})
        string(REGEX REPLACE "([0-9])" "0\\1" ver ${ver})
        set(${version} "0x${ver}")
    endif()
endmacro()
get_WIN32_WINNT(ver)
if(WIN32)
    add_definitions(-D_WIN32_WINNT=${ver} -DNOMINMAX)
    add_definitions(-DMSVC_EXPRESS -D_CRT_SECURE_NO_WARNINGS -DNGS_EXPORTS)
#     ADD_DEFINITIONS(-DNETGEN_ELTRANS) 
    
    # build convenience (aka object) libraries in windows)
    set(NGS_LIB_TYPE OBJECT)

    # omit openmp flags due to lacking support of openmp 3.0 in visual studio
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP /wd4067 /wd4081 /wd4018 /wd4373 /wd4592 /wd4275 /wd4522 /wd4910 /wd4244 /wd4800 /bigobj") #${OpenMP_CXX_FLAGS} 
else(WIN32)
    # build shared libraries in windows)
    set(NGS_LIB_TYPE SHARED)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
endif(WIN32)
if(APPLE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -undefined dynamic_lookup")
    add_definitions(-DMSG_NOSIGNAL=0)
#     set(CMAKE_SHARED_LIBRARY_SUFFIX ".so")
endif(APPLE)

#######################################################################
if (USE_GUI)
    find_package(TCL 8.5)
    if(TCL_FOUND)
        add_definitions(-DTCL)
        include_directories(${TCL_INCLUDE_PATH})
    else(TCL_FOUND)
        message( FATAL_ERROR "Tcl library not found!" )
    endif(TCL_FOUND)
endif (USE_GUI)

#######################################################################
if (USE_MPI)
    find_package(MPI REQUIRED)
    add_definitions(-DPARALLEL)
    include_directories(${MPI_CXX_INCLUDE_PATH})
    set(CMAKE_EXE_LINKER_FLAGS    "${CMAKE_EXE_LINKER_FLAGS}    ${MPI_CXX_LINK_FLAGS} -ldl")
    set(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} ${MPI_CXX_LINK_FLAGS} -ldl")
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${MPI_CXX_LINK_FLAGS} -ldl")
endif (USE_MPI)
# set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wl,--no-undefined")

#######################################################################
if (USE_LAPACK AND NOT USE_MKL)
    find_package(LAPACK)
    add_definitions(-DLAPACK)
endif (USE_LAPACK AND NOT USE_MKL)

#######################################################################
if (USE_VT)
    add_definitions(-DVTRACE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -t:inst manual")
endif (USE_VT)

#######################################################################
if (USE_CUDA)
    find_package(CUDA)
    set(CUDA_SEPARABLE_COMPILATION ON)
    set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS};-std=c++11;-O2;-DUSE_MYCOMPLEX;-DHD="__host__ __device__")
    set(CUDA_OPTIONS -Xcompiler -fPIC)
    SET (EXTRA_LIBS ${EXTRA_LIBS} cudart cublas cusparse)
endif (USE_CUDA)

#######################################################################
if (USE_MKL)
    set(MKL_MULTI_THREADED ON)
    set(MKL_STATIC OFF)
    set(MKL_SDL ON)
    set(USE_LAPACK ON)
    find_package(MKL REQUIRED)
    include_directories(${MKL_INCLUDE_DIRS})
    add_definitions(-DUSE_MKL -DUSE_PARDISO -DLAPACK)
    set( MKL_LIBRARIES "${MKL_MINIMAL_LIBRARIES}")
    if(NOT WIN32)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m64")
    endif(NOT WIN32)
endif (USE_MKL)

#######################################################################
if (USE_MUMPS)
    add_definitions(-DUSE_MUMPS)
    find_package(MUMPS REQUIRED)
    include_directories(${MUMPS_DIR}/include)
    link_directories(${MUMPS_DIR}/lib)
endif (USE_MUMPS)

#######################################################################
if (USE_PYTHON)
    set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH}" "${PROJECT_SOURCE_DIR}/cmake_modules/python")
    set(PYTHON_VERSION "3" CACHE STRING "Python version (only Python >= 3.0 supported)")
    if( PYTHON_VERSION VERSION_LESS 3 )
        message(FATAL_ERROR "NGSolve supports only Python 3")
    endif( PYTHON_VERSION VERSION_LESS 3 )
    find_package(PythonInterp ${PYTHON_VERSION} REQUIRED)
    find_package(PythonLibs ${PYTHON_VERSION}  REQUIRED)
    set(Boost_FIND_QUIETLY ON)
    if(NOT Boost_LIBRARIES)
        FIND_PACKAGE(Boost COMPONENTS python${PYTHON_VERSION_MAJOR})
    endif(NOT Boost_LIBRARIES)
    if(NOT Boost_LIBRARIES)
        find_package(Boost COMPONENTS python-py${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR})
    endif(NOT Boost_LIBRARIES)
    if(NOT Boost_LIBRARIES)
        find_package(Boost COMPONENTS python)
    endif(NOT Boost_LIBRARIES)
    if(NOT Boost_LIBRARIES)
        message(FATAL_ERROR "Could NOT find Boost-Python" )
    endif(NOT Boost_LIBRARIES)

    add_definitions(-DNGS_PYTHON)
    include_directories(${PYTHON_INCLUDE_DIRS})
    include_directories(${Boost_INCLUDE_DIRS})
    set(PYTHON_LIBS "${PYTHON_LIBRARIES};${Boost_LIBRARIES}")
    execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "from distutils.sysconfig import get_python_lib; print(get_python_lib(1,0,''))" OUTPUT_VARIABLE PYTHON_PACKAGES_INSTALL_DIR OUTPUT_STRIP_TRAILING_WHITESPACE)
    message("python install path: ${PYTHON_PACKAGES_INSTALL_DIR}")
endif (USE_PYTHON)

#######################################################################
if(USE_GOLD)
    add_definitions(-DGOLD -DJS)
    find_path(NGS_GOLD_DIR ngstd_gold.hpp)
    include_directories(${NGS_GOLD_DIR})
    if(NOT INTEL_MIC)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx2 -mfma")
    endif(NOT INTEL_MIC)
endif(USE_GOLD)

#######################################################################
if(INTEL_MIC)
    set(MKL_ARCH "mic")
    include(cmake_modules/mic.cmake)
else(INTEL_MIC)
    set(MKL_ARCH "intel64")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
endif(INTEL_MIC)

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

#######################################################################
# generate ngscxx compiler wrapper

get_directory_property( DirDefs COMPILE_DEFINITIONS )
set( NGSCXX_DEFS "" )
foreach( d ${DirDefs} )
    set( NGSCXX_DEFS "${NGSCXX_DEFS} -D${d}" )
endforeach()
message( STATUS "NGSCXX_DEFS: " ${NGSCXX_DEFS} )

get_directory_property(NGSCXX_INCLUDE_DIRS INCLUDE_DIRECTORIES)
foreach( d ${NGSCXX_INCLUDE_DIRS} )
    string(COMPARE EQUAL ${d} "${PROJECT_SOURCE_DIR}/include" IS_SOURCE_DIR)
    string(COMPARE EQUAL ${d} "${PROJECT_BINARY_DIR}" IS_BINARY_DIR)
    if(NOT IS_SOURCE_DIR AND NOT IS_BINARY_DIR)
        SET( NGSCXX_INCLUDE_FLAGS "${NGSCXX_INCLUDE_FLAGS} -I${d}" )
    endif(NOT IS_SOURCE_DIR AND NOT IS_BINARY_DIR)
endforeach()
list( REMOVE_DUPLICATES NGSCXX_INCLUDE_FLAGS )
message( STATUS "NGSCXX_INCLUDE_FLAGS: " ${NGSCXX_INCLUDE_FLAGS} )

if(NOT APPLE)
    set(ngscxx_extra_flags "${ngscxx_extra_flags} -Wl,--rpath=${INSTALL_DIR}/lib")  
endif(NOT APPLE)
string(TOUPPER "CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE}" name)
set(flags "${${name}} ${CMAKE_CXX_FLAGS}")
file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/ngscxx "${CMAKE_CXX_COMPILER} -I. ${flags} -L${INSTALL_DIR}/lib -fpic ${NGSCXX_DEFS} ${NGSCXX_INCLUDE_FLAGS} \$* ${ngscxx_extra_flags} \n")


#######################################################################
# CTest
enable_testing()
include(CTest)

#######################################################################
# Debian packager
set(CPACK_PACKAGE_VERSION ${PACKAGE_VERSION} )
set(CPACK_PACKAGING_INSTALL_PREFIX "${INSTALL_DIR}")
set(CPACK_PACKAGE_NAME ngsolve)

set(CPACK_SOURCE_GENERATOR "TGZ")
set(CPACK_SOURCE_IGNORE_FILES "/cmake/;/build/;/.gz/;~$;${CPACK_SOURCE_IGNORE_FILES}")

set(CPACK_GENERATOR "DEB")
execute_process(COMMAND dpkg --print-architecture OUTPUT_VARIABLE CPACK_DEBIAN_PACKAGE_ARCHITECTURE OUTPUT_STRIP_TRAILING_WHITESPACE)
set(CPACK_DEBIAN_PACKAGE_MAINTAINER Matthias Hochsteger <matthias.hochsteger@tuwien.ac.at>)
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Short description of ngsolve")
set(CPACK_PACKAGE_DESCRIPTION "Long description of ngsolve")
set(CPACK_DEBIAN_PACKAGE_DEPENDS "netgen (>= ${NGSOLVE_VERSION_MAJOR}.${NGSOLVE_VERSION_MINOR}), liblapacke")
# TODO (use dpkg -s <packagename> to retrieve version)

set(CPACK_DEBIAN_PACKAGE_SECTION Science)
set(CPACK_DEBIAN_PACKAGE_VERSION ${CPACK_PACKAGE_VERSION})
include(CPack)

#######################################################################
# uninstall target
configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake_modules/cmake_uninstall.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
    IMMEDIATE @ONLY)
add_custom_target(uninstall
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)

#######################################################################
# Configure message
# TODO: other message in case of failure
message("
------------------------------------------------------------------------
  ${PROJECT_NAME} ${NGSOLVE_VERSION}:  Automatic configuration OK.

  Build type: ${CMAKE_BUILD_TYPE}
  Flags: ${flags}

  Enabled functionality:

    GUI: ............... ${USE_GUI}
    LAPACK: ............ ${USE_LAPACK}
    MKL: ............... ${USE_MKL}
    PYTHON: ............ ${USE_PYTHON}
    MPI: ............... ${USE_MPI}
    CUDA: .............. ${USE_CUDA}
    MUMPS: ............. ${USE_MUMPS}
    GOLD: .............. ${USE_GOLD}
    INTEL_MIC: ......... ${INTEL_MIC}
    VTUNE: ............. ${USE_VTUNE}
    VAMPIRTRACE: ....... ${USE_VT}


  Building:

    Type 'make' to compile ${PROJECT_NAME}.

    Type 'make install' to install ${PROJECT_NAME}.

  Install directory:

    ${CMAKE_INSTALL_PREFIX}")
if(USE_MKL)
message("
  Please set the following environment variables:

    MKL_THREADING_LAYER=GNU")
endif(USE_MKL)
message("
------------------------------------------------------------------------
")


