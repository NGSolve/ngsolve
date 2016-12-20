if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "RelWithDebInfo" CACHE STRING INTERNAL)
endif(NOT CMAKE_BUILD_TYPE)

cmake_minimum_required(VERSION 2.8)
project(NGSolve)
set(NGSOLVE_VERSION_MAJOR 6)
set(NGSOLVE_VERSION_MINOR 2)
string(TIMESTAMP NGSOLVE_VERSION_PATCH "%y%U%w" )
# set(NGSOLVE_VERSION_PATCH 15440)
set(NGSOLVE_VERSION "${NGSOLVE_VERSION_MAJOR}.${NGSOLVE_VERSION_MINOR}-dev")
set(PACKAGE_VERSION "${NGSOLVE_VERSION_MAJOR}.${NGSOLVE_VERSION_MINOR}-${NGSOLVE_VERSION_PATCH}")
set(CPACK_PACKAGE_VERSION "${NGSOLVE_VERSION_MAJOR}.${NGSOLVE_VERSION_MINOR}-${NGSOLVE_VERSION_PATCH}")
set(CPACK_DEBIAN_PACKAGE_VERSION "${NGSOLVE_VERSION_MAJOR}.${NGSOLVE_VERSION_MINOR}.${NGSOLVE_VERSION_PATCH}")

if(NOT WIN32 AND NOT INTEL_MIC)
    option( USE_NATIVE_ARCH  "build which -march=native" ON)
endif(NOT WIN32 AND NOT INTEL_MIC)

option( USE_GUI          "build ngsolve with GUI" ON)
option( USE_PYTHON       "enable python shell" ON)
option( USE_LAPACK       "build ngsolve with LAPACK" ON)
option( USE_MPI          "enable mpi parallelization")
option( USE_VT           "enable vampirtrace")
option( USE_CUDA         "enable CUDA GPU support")
option( USE_MKL          "enable MKL")
option( USE_MUMPS        "enable sparse direct solver MUMPS")
option( USE_PARDISO      "enable pardiso sparse direct solver")
option( USE_UMFPACK      "enable umfpack sparse direct solver")
option( INTEL_MIC        "cross compile for intel xeon phi")
option( USE_VTUNE        "include vtune pause/resume numproc")
option( USE_NUMA         "compile with NUMA-aware code")
option( USE_CCACHE       "use ccache")
option( INSTALL_DEPENDENCIES "install dependencies like netgen or solver libs, useful for packaging" OFF )

set(NETGEN_INTERNAL_INCLUDE_DIR NOTFOUND CACHE INTERNAL "Netgen include directories")
set(BUILD_NETGEN OFF CACHE INTERNAL "is set if netgen is built as external project")
#######################################################################
if(NETGEN_SOURCE_DIR)
    find_path(NETGEN_INCLUDE_DIR nginterface_v2_impl.hpp PATHS ${NETGEN_SOURCE_DIR} ${NETGEN_SOURCE_DIR}/libsrc/include)
    if (NOT NETGEN_INCLUDE_DIR)
      message(FATAL_ERROR "Could not find Netgen source files in NETGEN_SOURCE_DIR, which was set to ${NETGEN_SOURCE_DIR}")
    endif()
    set(NETGEN_INTERNAL_INCLUDE_DIR
      "${NETGEN_INCLUDE_DIR}"
      "${NETGEN_INCLUDE_DIR}/../general"
      "${NETGEN_INCLUDE_DIR}/../visualization"
    )
else(NETGEN_SOURCE_DIR)
  message(STATUS "Use Netgen from submodule, updating modules...")
  execute_process(COMMAND git submodule update --init --recursive WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
  add_custom_target(check_submodules_start ALL cmake -P cmake_modules/check_submodules.cmake WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} )
  add_custom_target(check_submodules_stop ALL cmake -P cmake_modules/check_submodules.cmake WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} DEPENDS ngs ngslib netgen ngpy occvis)
  set(BUILD_NETGEN ON)
  set(NETGEN_INTERNAL_INCLUDE_DIR
    "${CMAKE_CURRENT_SOURCE_DIR}/external_dependencies/netgen/libsrc/include"
    "${CMAKE_CURRENT_SOURCE_DIR}/external_dependencies/netgen/libsrc/general"
    "${CMAKE_CURRENT_SOURCE_DIR}/external_dependencies/netgen/libsrc/visualization"
  )
endif(NETGEN_SOURCE_DIR)

#######################################################################
include (CheckIncludeFiles)
check_include_files (dlfcn.h HAVE_DLFCN_H)
include (CheckCXXSymbolExists)
# set(CMAKE_REQUIRED_INCLUDES cxxabi.h)
check_cxx_symbol_exists ( abi::__cxa_demangle "cxxabi.h" HAVE_CXA_DEMANGLE )
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/config.h.cmake ${CMAKE_CURRENT_BINARY_DIR}/config.h)

install (FILES ${CMAKE_CURRENT_BINARY_DIR}/config.h DESTINATION include COMPONENT ngsolve_devel )


add_subdirectory(external_dependencies)

include_directories(AFTER ${NETGEN_INTERNAL_INCLUDE_DIR})
#######################################################################

add_definitions(-DHAVE_NETGEN_SOURCES)
#######################################################################
if(USE_CCACHE)
  find_program(CCACHE_FOUND ccache)
  if(CCACHE_FOUND)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE ccache)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK ccache)
  endif(CCACHE_FOUND)
endif(USE_CCACHE)

#######################################################################
if(INTEL_MIC)
    set(MKL_ARCH "mic")
    include(cmake_modules/mic.cmake)
endif(INTEL_MIC)

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
if(APPLE)
  set(CMAKE_INSTALL_PREFIX "${INSTALL_DIR}/Netgen.app/Contents/Resources" CACHE INTERNAL "Prefix prepended to install directories" FORCE)
else(APPLE)
  set(CMAKE_INSTALL_PREFIX "${INSTALL_DIR}" CACHE INTERNAL "Prefix prepended to install directories" FORCE) 
endif(APPLE)

# add_library(interface SHARED IMPORTED)
# set_property(TARGET interface PROPERTY IMPORTED_LOCATION "${INSTALL_DIR}/lib/libinterface.so")
#######################################################################
# append install paths of software in non-standard paths (e.g. openmpi, metis, intel mkl, ...)
# cmake -DUSE_MPI=ON -DCMAKE_PREFIX_PATH="/opt/openmpi165;/opt/metis51" ../
set(ADDITIONAL_PATHS "" CACHE PATH "List of paths to additional libraries in non-standard locations, seperated by ';'")
if (ADDITIONAL_PATHS)
  set(CMAKE_PREFIX_PATH ${ADDITIONAL_PATHS})
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

include_directories ("${PROJECT_SOURCE_DIR}/include")
include_directories ("${PROJECT_BINARY_DIR}")
# include_directories (AFTER "${CMAKE_INSTALL_PREFIX}/include")

link_directories ("${CMAKE_INSTALL_PREFIX}/${ngs_install_dir_lib}")

install (FILES ngsolve.tcl DESTINATION ${ngs_install_dir_bin} COMPONENT ngsolve )

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DDEBUG -Wall")

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

    # build convenience (aka object) libraries in windows)
    set(NGS_LIB_TYPE OBJECT)

    if(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
        # Intel Compiler
        set(CMAKE_CXX_FLAGS "/Od ${CMAKE_CXX_FLAGS} /MP /Qstd=c++14 /Qopt-report-file:report.log")
    else()
        # MS Visual Compiler
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP /w /bigobj")
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
    if (USE_PYTHON)
      find_package(MPI4PY REQUIRED)
      include_directories(${MPI4PY_INCLUDE_DIR})
    endif (USE_PYTHON)
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
    set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS};-std=c++14;-O2;-DUSE_MYCOMPLEX;-DHD="__host__ __device__")
    set(CUDA_OPTIONS -Xcompiler -fPIC)
    SET (EXTRA_LIBS ${EXTRA_LIBS} cudart cublas cusparse)
endif (USE_CUDA)

#######################################################################
if (USE_MKL)
    set(MKL_MULTI_THREADED ON)

    set(MKL_STATIC OFF CACHE BOOL "Link static MKL")
    set(MKL_SDL ON CACHE BOOL "Link single dynamic MKL lib")

    set(USE_LAPACK ON)
    find_package(MKL REQUIRED)
    include_directories(${MKL_INCLUDE_DIRS})
    add_definitions(-DUSE_MKL -DUSE_PARDISO -DLAPACK)

    if(USE_MUMPS)
        # include scalapack
        set( LAPACK_LIBRARIES "${MKL_LIBRARIES}")
    else(USE_MUMPS)
        set( LAPACK_LIBRARIES "${MKL_MINIMAL_LIBRARIES}")
    endif(USE_MUMPS)

    if("${MKL_INTERFACE_LAYER}" STREQUAL "_ilp64")
        add_definitions(-DMKL_ILP64)
    endif()
endif (USE_MKL)

#######################################################################
if (USE_MUMPS)
    add_definitions(-DUSE_MUMPS)
    find_package(MUMPS REQUIRED)
    set(MUMPS_INCLUDE_DIR ${MUMPS_DIR}/include)
    include_directories(${MUMPS_INCLUDE_DIR})
    link_directories(${MUMPS_DIR}/lib)
endif (USE_MUMPS)

#######################################################################
if (USE_PYTHON)
  find_path(PYBIND_INCLUDE_DIR pybind11/pybind11.h ${NETGEN_INTERNAL_INCLUDE_DIR}/../../external_dependencies/pybind11/include)
    if( NOT PYBIND_INCLUDE_DIR )
      message(FATAL_ERROR "Could NOT find pybind11!")
    endif( NOT PYBIND_INCLUDE_DIR )
    message("-- Found Pybind11: ${PYBIND_INCLUDE_DIR}")
    include_directories(BEFORE ${PYBIND_INCLUDE_DIR})
    install(DIRECTORY ${PYBIND_INCLUDE_DIR}/pybind11 DESTINATION include COMPONENT ngsolve_devel)
    install(FILES ${PYBIND_INCLUDE_DIR}/../LICENSE DESTINATION include/pybind11 COMPONENT ngsolve_devel)
    set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH}" "${PROJECT_SOURCE_DIR}/cmake_modules/python")
    set(PYTHON_VERSION "3" CACHE STRING "Python version (only Python >= 3.0 supported)")
    set(Python_ADDITIONAL_VERSIONS 3.5)
    if( PYTHON_VERSION VERSION_LESS 3 )
        message(FATAL_ERROR "NGSolve supports only Python 3")
    endif( PYTHON_VERSION VERSION_LESS 3 )
    find_package(PythonInterp ${PYTHON_VERSION} REQUIRED)
    find_package(PythonLibs ${PYTHON_VERSION}  REQUIRED)
    add_definitions(-DNGS_PYTHON)

    include_directories(${PYTHON_INCLUDE_DIRS})
    set(PYTHON_LIBS "${PYTHON_LIBRARIES}")
    execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "from distutils.sysconfig import get_python_lib; print(get_python_lib(1,0,''))" OUTPUT_VARIABLE PYTHON_PACKAGES_INSTALL_DIR OUTPUT_STRIP_TRAILING_WHITESPACE)
endif (USE_PYTHON)

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
  add_definitions(-DUSE_UMFPACK)

  if(NOT UMFPACK_INCLUDE_DIR)
    find_path(UMFPACK_INCLUDE_DIR umfpack.h PATH_SUFFIXES suitesparse)
  endif()

  if(NOT UMFPACK_LIBRARIES)
    find_library(UMFPACK_LIB NAMES libumfpack.a )
    if(UMFPACK_LIB)
      find_library(UMFPACK_AMD_LIB NAMES libamd.a)
      find_library(UMFPACK_COLAMD_LIB NAMES libcolamd.a)
      find_library(UMFPACK_CONFIG_LIB NAMES libsuitesparseconfig.a)
      find_library(UMFPACK_CHOLMOD_LIB NAMES libcholmod.a)
    else(UMFPACK_LIB)
      find_library(UMFPACK_LIB NAMES libumfpack.so)
    endif(UMFPACK_LIB)
  endif(NOT UMFPACK_LIBRARIES)

  set(UMFPACK_LIBRARIES ${UMFPACK_LIB} ${UMFPACK_AMD_LIB} ${UMFPACK_COLAMD_LIB} ${UMFPACK_CHOLMOD_LIB} ${UMFPACK_CONFIG_LIB} CACHE STRING "List of UMFPACK libraries")

  if(NOT UMFPACK_LIBRARIES OR NOT UMFPACK_INCLUDE_DIR)
    message(FATAL_ERROR "Could NOT find UMFPACK")
  else()
    message(STATUS "Found Umfpack include:   ${UMFPACK_INCLUDE_DIR}")
    message(STATUS "Found Umfpack libraries: ${UMFPACK_LIBRARIES}")
  endif()

  include_directories(${UMFPACK_INCLUDE_DIR})

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

if(BUILD_NETGEN)
  add_dependencies(ngstd check_submodules_start)
  add_dependencies(ngmg check_submodules_start)
endif(BUILD_NETGEN)

#######################################################################
# generate ngscxx compiler wrapper

# compile flags and preprocessor definitions
get_directory_property( ngscxx_defines COMPILE_DEFINITIONS )
string(TOUPPER "CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE}" name)

foreach( d ${ngscxx_defines} )
  if(NOT ${d} MATCHES "HAVE_NETGEN_SOURCES")
    if(WIN32)
      if(NOT ${d} STREQUAL "NGS_EXPORTS")
        set( ngscxx_define_flags "${ngscxx_define_flags} /D${d}" )
      endif()
    else()
      set( ngscxx_define_flags "${ngscxx_define_flags} -D${d}" )
    endif()
  endif()
endforeach()

set(ngscxx_flags "${${name}} ${CMAKE_CXX_FLAGS} -fPIC ${ngscxx_define_flags}")

# include directories
set(ngscxx_include_dirs ${INSTALL_DIR}/include ${TCL_INCLUDE_PATH} ${MPI_CXX_INCLUDE_PATH} ${MKL_INCLUDE_DIRS} ${MUMPS_INCLUDE_DIR} ${PYTHON_INCLUDE_DIRS} ${UMFPACK_INCLUDE_DIR} ${VTUNE_INCLUDE_DIR})
if(NOT INSTALL_DEPENDENCIES) 
    set(ngscxx_include_dirs ${ngscxx_include_dirs} ${NETGEN_INTERNAL_INCLUDE_DIR})
endif(NOT INSTALL_DEPENDENCIES) 


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



#######################################################################
# CTest
enable_testing()
include(CTest)
add_custom_target(performance_tests ctest -L performance -D Experimental)

#######################################################################
# Debian packager
set(CPACK_PACKAGE_VERSION ${PACKAGE_VERSION} )
set(CPACK_PACKAGE_NAME ngsolve)

set(CPACK_SOURCE_GENERATOR "TGZ")
set(CPACK_SOURCE_IGNORE_FILES "/cmake/;/build/;/.gz/;~$;${CPACK_SOURCE_IGNORE_FILES}")

set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "NGSolve Finite Element Library")
set(CPACK_PACKAGE_DESCRIPTION "NGSolve is a general purpose Finite Element Library on top of Netgen. With the basic library one can solve heat flow equations, Maxwell equations, and solid mechanical problems. Several add-ons are available for particular application classes.") 
if(UNIX)
    set(CPACK_PACKAGING_INSTALL_PREFIX "/opt/netgen")
    set(CPACK_GENERATOR "DEB")
    execute_process(COMMAND dpkg --print-architecture OUTPUT_VARIABLE CPACK_DEBIAN_PACKAGE_ARCHITECTURE OUTPUT_STRIP_TRAILING_WHITESPACE)
    set(CPACK_DEBIAN_PACKAGE_MAINTAINER "Matthias Hochsteger <matthias.hochsteger@tuwien.ac.at>")
    set(CPACK_DEBIAN_PACKAGE_DEPENDS "netgen (>= ${PACKAGE_VERSION}), liblapacke")
# TODO (use dpkg -s <packagename> to retrieve version)

    set(CPACK_DEBIAN_PACKAGE_SECTION Science)
    set(CPACK_PACKAGE_NAME ${CPACK_PACKAGE_NAME})
    set(CPACK_DEBIAN_PACKAGE_NAME ${CPACK_PACKAGE_NAME})
    set(CPACK_PACKAGE_FILE_NAME "ngsolve-${PACKAGE_VERSION}_${CPACK_DEBIAN_PACKAGE_ARCHITECTURE}")
endif(UNIX)
if(WIN32)
    set(CPACK_GENERATOR "NSIS")
    set(CPACK_PACKAGE_INSTALL_DIRECTORY "netgen 6.1-dev")
    set(CPACK_INSTALL_PREFIX "")
    set(CPACK_PACKAGING_INSTALL_PREFIX "")
endif(WIN32)

if(APPLE)
    # make sure to package netgen with all its dependencies
    if(INSTALL_DEPENDENCIES)
        file(STRINGS ${CMAKE_CURRENT_BINARY_DIR}/../netgen/install_manifest.txt netgen_install_files)
        foreach( f ${netgen_install_files} )
          string(REGEX REPLACE "^.*/(lib|share|doc|include|bin|../Frameworks)/((([^/])*/)*)[^/]*" "\\1/\\2" fout ${f})
          if(NOT IS_DIRECTORY ${f})
            if(${fout} STREQUAL "bin/")
              install(PROGRAMS ${f} DESTINATION ${fout})
            else()
              install(FILES ${f} DESTINATION ${fout})
            endif()
          endif()
        endforeach()
    endif(INSTALL_DEPENDENCIES)
#     set(CPACK_GENERATOR "Bundle")
#     set(CPACK_BUNDLE_NAME "NGSuite")
#     set(CPACK_PACKAGE_FILE_NAME "ngsuite-${PACKAGE_VERSION}")
#     set(CPACK_BUNDLE_ICON ${CMAKE_CURRENT_SOURCE_DIR}/logo.icns)
    #   set(CPACK_PACKAGE_ICON ${CMAKE_CURRENT_SOURCE_DIR}/logo.icns)            
    install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/logo.icns DESTINATION . RENAME Netgen.icns)


    # create some auxiliary files
    set(mac_startup ${CMAKE_CURRENT_BINARY_DIR}/startup.sh)
    file(WRITE ${mac_startup} "#!/bin/sh\n")
    file(APPEND ${mac_startup} "NGSuite_BUNDLE=\"`echo \"$0\" | sed -e 's/\\/Contents\\/Resources\\/bin\\/startup.sh//'`\"\n")
    file(APPEND ${mac_startup} "NGSuite_RESOURCES=\"$NGSuite_BUNDLE/Contents/Resources\"\n")
    file(APPEND ${mac_startup} "export DYLD_LIBRARY_PATH=$NGSuite_RESOURCES/lib:$DYLD_LIBRARY_PATH\n")
    file(APPEND ${mac_startup} "export PYTHONPATH=$NGSuite_RESOURCES/${PYTHON_PACKAGES_INSTALL_DIR}:$PYTHONPATH")
    file(APPEND ${mac_startup} "export NETGENDIR=$NGSuite_RESOURCES/bin\n")
    file(APPEND ${mac_startup} "cd $NGSuite_RESOURCES\n")
    file(APPEND ${mac_startup} "$NGSuite_RESOURCES/bin/netgen\n")
#     install(PROGRAMS ${mac_startup} DESTINATION bin)

    set(mac_ngsuite ${CMAKE_CURRENT_BINARY_DIR}/ngsuite.sh)
    file(WRITE ${mac_ngsuite} "#!/bin/sh\n")
    file(APPEND ${mac_ngsuite} "NGSuite_BUNDLE=\"`echo \"$0\" | sed -e 's/\\/Contents\\/MacOS\\/NGSuite//'`\"\n")
    file(APPEND ${mac_ngsuite} "NGSuite_RESOURCES=\"$NGSuite_BUNDLE/Contents/Resources\"\n")
    file(APPEND ${mac_ngsuite} "open -a /Applications/Utilities/Terminal.app $NGSuite_RESOURCES/bin/startup.sh\n")
#     set(CPACK_BUNDLE_STARTUP_COMMAND ${mac_ngsuite})

    set(mac_plist ${CMAKE_CURRENT_BINARY_DIR}/Info.plist)
    file(WRITE ${mac_plist} "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n")
    file(APPEND ${mac_plist} "<!DOCTYPE plist PUBLIC \"-//Apple Computer//DTD PLIST 1.0//EN\" \"http://www.apple.com/DTDs/PropertyList-1.0.dtd\">\n")
    file(APPEND ${mac_plist} "<plist version=\"1.0\">\n")
    file(APPEND ${mac_plist} "<dict>\n")
    file(APPEND ${mac_plist} "	<key>CFBundleDevelopmentRegion</key>\n")
    file(APPEND ${mac_plist} "	<string>English</string>\n")
    file(APPEND ${mac_plist} "	<key>CFBundleExecutable</key>\n")
    file(APPEND ${mac_plist} "	<string>NGSuite</string>\n")
    file(APPEND ${mac_plist} "	<key>CFBundleIconFile</key>\n")
    file(APPEND ${mac_plist} "	<string>NGSuite.icns</string>\n")
    file(APPEND ${mac_plist} "</dict>\n")
    file(APPEND ${mac_plist} "</plist>\n")
#     set(CPACK_BUNDLE_PLIST ${mac_plist})
#    set(CPACK_DMG_FORMAT "UDRW")
#    set(CPACK_DMG_DS_STORE .DS_Store_dmg)
#    set(CPACK_DMG_BACKGROUND_IMAGE sculpture312.png)
    # read install manifest of netgen
endif(APPLE)

include(CPack)

#######################################################################
# uninstall target
configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake_modules/cmake_uninstall.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
    IMMEDIATE @ONLY)
# add_custom_target(uninstall
#     COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cmake_modules/NGSolveConfig.cmake.in
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

    GUI: ............... ${USE_GUI}
    LAPACK: ............ ${USE_LAPACK}
    MKL: ............... ${USE_MKL}
    PYTHON: ............ ${USE_PYTHON}
    MPI: ............... ${USE_MPI}
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


