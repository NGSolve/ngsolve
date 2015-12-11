if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "RelWithDebInfo" CACHE STRING INTERNAL)
endif(NOT CMAKE_BUILD_TYPE)

cmake_minimum_required(VERSION 2.8.9)
project(Netgen)
set(NETGEN_VERSION_MAJOR 6)
set(NETGEN_VERSION_MINOR 1)
string(TIMESTAMP NETGEN_VERSION_PATCH "%y%U%w" )
set(NETGEN_VERSION "${NETGEN_VERSION_MAJOR}.${NETGEN_VERSION_MINOR}-dev")
set(PACKAGE_VERSION "${NETGEN_VERSION_MAJOR}.${NETGEN_VERSION_MINOR}-${NETGEN_VERSION_PATCH}")
set(CPACK_PACKAGE_VERSION "${PACKAGE_VERSION}")
# set(CPACK_DEBIAN_PACKAGE_VERSION "${PACKAGE_VERSION}"

if(NOT WIN32)
    option( USE_NATIVE_ARCH  "build which -march=native" ON)
endif(NOT WIN32)

option( USE_GUI     "don't build netgen with GUI" ON )
option( USE_PYTHON  "link boost python, necessary for ngsolve" ON )
option( USE_MPI     "enable mpi parallelization" OFF )
option( USE_OCC     "(not supported) compile with OpenCascade geometry kernel" OFF)
option( USE_JPEG    "enable snapshots using library libjpeg" OFF )
option( USE_MPEG    "enable video recording with FFmpeg, uses libavcodec" OFF )
option( INTEL_MIC        "cross compile for intel xeon phi")
option( INSTALL_PROFILES "install environment variable settings to /etc/profile.d" OFF )
option( INSTALL_DEPENDENCIES "install dependencies like tcl, useful for packaging" OFF )
option( USE_CCACHE       "use ccache")


#######################################################################
if(USE_CCACHE)
  find_program(CCACHE_FOUND ccache)
  if(CCACHE_FOUND)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE ccache)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK ccache)
  endif(CCACHE_FOUND)
endif(USE_CCACHE)

#######################################################################
if(USE_NATIVE_ARCH)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
endif(USE_NATIVE_ARCH)

#######################################################################
if(INTEL_MIC)
    set(MKL_ARCH "mic")
    include(cmake_modules/mic.cmake)
else(INTEL_MIC)
    set(MKL_ARCH "intel64")
endif(INTEL_MIC)

#######################################################################
# Append install paths of software in non-standard paths (e.g. openmpi, metis, intel mkl, ...)
# cmake -DUSE_MPI=ON -DADDITIONAL_PATHS="/opt/openmpi165;/opt/metis51" ../
set(ADDITIONAL_PATHS "" CACHE PATH "List of paths to additional libraries in non-standard locations, seperated by ';'")
if (ADDITIONAL_PATHS)
    list(REVERSE CMAKE_SYSTEM_PREFIX_PATH)
    list(APPEND CMAKE_SYSTEM_PREFIX_PATH ${ADDITIONAL_PATHS})
    list(REVERSE CMAKE_SYSTEM_PREFIX_PATH)
endif (ADDITIONAL_PATHS)

#######################################################################
# install dir
if(NOT INSTALL_DIR)
    set(INSTALL_DIR "/opt/netgen" CACHE PATH "Prefix prepended to install directories")
endif(NOT INSTALL_DIR)

set(CMAKE_INSTALL_PREFIX "${INSTALL_DIR}" CACHE INTERNAL "Prefix prepended to install directories" FORCE)

#######################################################################
# use rpath
SET(CMAKE_SKIP_BUILD_RPATH  FALSE)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE) 
SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
IF("${isSystemDir}" STREQUAL "-1")
   SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
ENDIF("${isSystemDir}" STREQUAL "-1")

#######################################################################
# build options
add_definitions(-DUSE_TIMEOFDAY)

include_directories ("${PROJECT_SOURCE_DIR}/include")
include_directories ("${PROJECT_SOURCE_DIR}/libsrc/include")
include_directories ("${PROJECT_BINARY_DIR}")

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH}" "${PROJECT_SOURCE_DIR}/cmake_modules")
set(ng_install_dir RUNTIME DESTINATION bin COMPONENT netgen LIBRARY DESTINATION lib COMPONENT netgen_devel ARCHIVE DESTINATION lib COMPONENT netgen_devel)

include (CheckIncludeFiles)
check_include_files (dlfcn.h HAVE_DLFCN_H)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/config.h.cmake ${CMAKE_CURRENT_BINARY_DIR}/config.h)
add_definitions(-DHAVE_CONFIG_H)

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

macro(get_dll_from_lib dll_path lib_path)
        get_filename_component(parent_lib_path ${lib} DIRECTORY)
        get_filename_component(lib_name ${lib} name)
endmacro()

if(WIN32)
    get_WIN32_WINNT(ver)
    add_definitions(-D_WIN32_WINNT=${ver})
    set(CMAKE_MFC_FLAG 0)
#     add_definitions(-DNGINTERFACE_EXPORTS)
#     add_definitions(-DNGLIB_EXPORTS)
    add_definitions(-DMSVC_EXPRESS -D_CRT_SECURE_NO_WARNINGS -DHAVE_STRUCT_TIMESPEC)
    # build convenience (aka object) libraries in windows)
    set(NG_LIB_TYPE OBJECT)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP /wd4244 /wd4800")
else(WIN32)
    # build shared libraries
    set(NG_LIB_TYPE SHARED)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
endif(WIN32)
if(APPLE)
#    set(MACOSX_BUNDLE ON)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -undefined dynamic_lookup")
#     set(CMAKE_SHARED_LIBRARY_SUFFIX ".so")
endif(APPLE)

#######################################################################
if(NOT INTEL_MIC AND NOT WIN32)
    find_package(OpenMP)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif(NOT INTEL_MIC AND NOT WIN32)

#######################################################################
find_package(ZLIB REQUIRED)
include_directories(${ZLIB_INCLUDE_DIRS})

#######################################################################
if (USE_GUI)
    set(CMAKE_THREAD_PREFER_PTHREAD ON)
    find_package(TclStub 8.5 REQUIRED)
    find_package(Threads REQUIRED)
    find_package(X11 REQUIRED)
    if(APPLE)
        set(OPENGL_glu_LIBRARY "/opt/X11/lib/libGLU.dylib") 
        set(OPENGL_gl_LIBRARY "/opt/X11/lib/libGL.dylib") 
    endif(APPLE)
    find_package(OpenGL REQUIRED)
    find_library(LIBPTHREAD NAMES pthread HINTS ${MY_LIB_DIR} )

    add_definitions(-DTCL -DOPENGL)
    include_directories(${TCL_INCLUDE_PATH})
#     include_directories(${TK_INCLUDE_PATH})

    if(WIN32)
        get_filename_component(MY_LIB_DIR ${TK_LIBRARY} DIRECTORY)
        find_library(LIBTOGL NAMES togl togl17 HINTS ${MY_LIB_DIR} )
        install( DIRECTORY "${MY_LIB_DIR}/tcl8.5"   DESTINATION lib COMPONENT netgen )
        install( DIRECTORY "${MY_LIB_DIR}/tk8.5"    DESTINATION lib COMPONENT netgen )
        install( DIRECTORY "${MY_LIB_DIR}/tix8.4.3" DESTINATION lib COMPONENT netgen )
        install( DIRECTORY "${MY_LIB_DIR}/../bin" DESTINATION . COMPONENT netgen )
    else(WIN32)
        set(LIBTOGL togl)
        ADD_DEFINITIONS(-DTOGL_X11)
    endif(WIN32)
endif (USE_GUI)

#######################################################################
if (USE_PYTHON)
    set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH}" "${PROJECT_SOURCE_DIR}/cmake_modules/python")
    set(PYTHON_VERSION "3" CACHE STRING "")
    set(Python_ADDITIONAL_VERSIONS 3.5)
    if( PYTHON_VERSION VERSION_LESS 3 )
        message(FATAL_ERROR "NGSolve supports only Python 3")
    endif( PYTHON_VERSION VERSION_LESS 3 )
    if( PYTHON_VERSION STREQUAL "3")
        find_package(PythonInterp ${PYTHON_VERSION} REQUIRED)
        find_package(PythonLibs ${PYTHON_VERSION}  REQUIRED)
    else( PYTHON_VERSION STREQUAL "3")
        find_package(PythonInterp ${PYTHON_VERSION} EXACT REQUIRED)
        find_package(PythonLibs ${PYTHON_VERSION}  EXACT REQUIRED)
    endif( PYTHON_VERSION STREQUAL "3")
    set(Boost_FIND_QUIETLY ON)
    if(NOT Boost_LIBRARIES)
        find_package(Boost COMPONENTS python${PYTHON_VERSION_MAJOR})
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

    add_definitions(-DNG_PYTHON)

    if(WIN32 AND ${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel")
        # fix problems with boosts config/auto_link.hpp
        add_definitions(-DBOOST_LIB_TOOLSET="iw")
    endif()

    include_directories(${PYTHON_INCLUDE_DIRS})
    include_directories(${Boost_INCLUDE_DIRS})
    set(PYTHON_LIBS "${PYTHON_LIBRARIES};${Boost_LIBRARIES}")
    execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "from distutils.sysconfig import get_python_lib; print(get_python_lib(1,0,''))" OUTPUT_VARIABLE PYTHON_PACKAGES_INSTALL_DIR OUTPUT_STRIP_TRAILING_WHITESPACE)
    message("python install path: ${PYTHON_PACKAGES_INSTALL_DIR}")
endif (USE_PYTHON)

#######################################################################
if (USE_MPI)
    find_package(MPI REQUIRED)
    find_package(METIS REQUIRED)
    add_definitions(-DPARALLEL -DMETIS)
    include_directories(${MPI_CXX_INCLUDE_PATH})
    include_directories(${METIS_INCLUDE_DIR})
endif (USE_MPI)

#######################################################################
if (USE_OCC)
    find_package(OpenCasCade REQUIRED)
    add_definitions(-DOCCGEOMETRY -D_OCC64)

	add_definitions(-DHAVE_IOSTREAM -DHAVE_IOSTREAM_H -DHAVE_LIMITS -DHAVE_LIMITS_H -DHAVE_IOMANIP -DHAVE_IOMANIP_H)
    include_directories(${OCC_INCLUDE_DIR})
endif (USE_OCC)

#######################################################################
if (USE_JPEG)
    find_package(JPEG REQUIRED)
    add_definitions(-DJPEGLIB)
    include_directories(${JPEG_INCLUDE_DIR})
endif (USE_JPEG)

#######################################################################
if (USE_MPEG)
    find_package(FFMPEG REQUIRED)
    add_definitions(-DFFMPEG -D__STDC_CONSTANT_MACROS)
    include_directories(${FFMPEG_INCLUDE_DIR})
endif (USE_MPEG)

#######################################################################
if(INSTALL_PROFILES)
    file(WRITE  ${CMAKE_CURRENT_BINARY_DIR}/netgen.sh "#!/bin/sh\n")
    file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/netgen.sh "export PATH=${INSTALL_DIR}/bin:$PATH\n")
    file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/netgen.sh "export NETGENDIR=${INSTALL_DIR}/bin\n")
    file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/netgen.sh "export PYTHONPATH=${INSTALL_DIR}/${PYTHON_PACKAGES_INSTALL_DIR}:.\n")
    file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/netgen.sh "export LD_LIBRARY_PATH=${INSTALL_DIR}/lib:.\n")
    install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/netgen.sh DESTINATION /etc/profile.d)

    string(ASCII 27 Esc)
    file(WRITE  ${CMAKE_CURRENT_BINARY_DIR}/postinst "#!/bin/sh\n")
    file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/postinst "echo \"${Esc}[0;31mPlease log out and in again or do 'source /etc/profile.d/netgen.sh' to load the correct environment variables!${Esc}[m\"")
    set( CPACK_DEBIAN_PACKAGE_CONTROL_EXTRA "${CMAKE_CURRENT_BINARY_DIR}/postinst")
endif(INSTALL_PROFILES)

#######################################################################
if(INSTALL_DEPENDENCIES)
  include (InstallRequiredSystemLibraries)
  message("***********************************************")
  find_library(LIBGOMP NAMES gomp.1)
  find_library(LIBSTDCXX NAMES stdc++.6)
  message("gomp: ${LIBGOMP}")
  message("stdc++: ${LIBSTDCXX}")
  get_filename_component(MY_LIB_DIR ${TK_LIBRARY} DIRECTORY)
  message("my_tcl lib dir: ${MY_LIB_DIR}")
  install( DIRECTORY "${MY_LIB_DIR}"   DESTINATION . COMPONENT netgen )
  install( FILES "${Boost_LIBRARIES}" ${LIBGOMP} ${LIBSTDCXX} DESTINATION lib COMPONENT netgen )
endif(INSTALL_DEPENDENCIES)
#######################################################################

add_subdirectory(libsrc)
add_subdirectory(ng)
add_subdirectory(tutorials)
if (USE_PYTHON)
  add_subdirectory(python)
endif (USE_PYTHON)
add_subdirectory(py_tutorials)
add_subdirectory(doc)
add_subdirectory(windows)
add_subdirectory(nglib)


#######################################################################
# CTest
enable_testing()
include(CTest)

#######################################################################
# Debian packager

set(CPACK_SOURCE_GENERATOR "TGZ")
set(CPACK_SOURCE_IGNORE_FILES "/cmake/;/build/;/.gz/;~$;${CPACK_SOURCE_IGNORE_FILES}")
set(CPACK_PACKAGE_VERSION ${PACKAGE_VERSION} )

set(CPACK_PACKAGE_NAME netgen)
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "automatic 3d tetrahedral mesh generator")
set(CPACK_PACKAGE_DESCRIPTION "NETGEN is an automatic 3d tetrahedral mesh generator. It accepts input from constructive solid geometry (CSG) or boundary representation (BRep) from STL file format. The connection to a geometry kernel allows the handling of IGES and STEP files. NETGEN contains modules for mesh optimization and hierarchical mesh refinement. Netgen is open source based on the LGPL license. It is available for Unix/Linux and Windows.")

if(UNIX)
    set(CPACK_PACKAGING_INSTALL_PREFIX "/opt/netgen")
    execute_process(COMMAND grep CODENAME /etc/lsb-release OUTPUT_VARIABLE temp OUTPUT_STRIP_TRAILING_WHITESPACE)
    execute_process(COMMAND dpkg --print-architecture OUTPUT_VARIABLE CPACK_DEBIAN_PACKAGE_ARCHITECTURE OUTPUT_STRIP_TRAILING_WHITESPACE)
    if(temp)
        set(CPACK_GENERATOR "DEB")
        string(SUBSTRING ${temp} 17 -1 UBUNTU_VERSION)
        message("ubuntu version: ${UBUNTU_VERSION}")

        set(BOOST_PACKAGE "libboost-python${Boost_MAJOR_VERSION}.${Boost_MINOR_VERSION}.${Boost_SUBMINOR_VERSION}")
        set(CPACK_DEBIAN_PACKAGE_DEPENDS "python3, ${BOOST_PACKAGE}, libtk8.5, libtcl8.5, tix, libxmu6")
        execute_process(COMMAND dpkg --print-architecture OUTPUT_VARIABLE CPACK_DEBIAN_PACKAGE_ARCHITECTURE OUTPUT_STRIP_TRAILING_WHITESPACE)
        set(CPACK_DEBIAN_PACKAGE_MAINTAINER "Matthias Hochsteger <matthias.hochsteger@tuwien.ac.at>")
        if(USE_MPI)
            set(CPACK_DEBIAN_PACKAGE_DEPENDS "${CPACK_DEBIAN_PACKAGE_DEPENDS}, libmetis5, openmpi-bin")
            set(CPACK_PACKAGE_NAME "${CPACK_PACKAGE_NAME}_mpi")
        endif(USE_MPI)
        set(CPACK_DEBIAN_PACKAGE_SECTION Science)
        set(CPACK_DEBIAN_PACKAGE_NAME ${CPACK_PACKAGE_NAME})
        set(CPACK_PACKAGE_FILE_NAME "netgen-${PACKAGE_VERSION}_${CPACK_DEBIAN_PACKAGE_ARCHITECTURE}")   
    endif(temp)
endif(UNIX)

if(APPLE)
    set(CPACK_GENERATOR "Bundle")
    set(CPACK_BUNDLE_NAME "NGSuite")
    set(CPACK_BUNDLE_ICON ${CMAKE_CURRENT_SOURCE_DIR}/../netgen/windows/netgen.ico)

    # create some auxiliary files
    set(mac_startup ${CMAKE_CURRENT_BINARY_DIR}/startup.sh)
    file(WRITE ${mac_startup} "#!/bin/sh\n")
    file(APPEND ${mac_startup} "NGSuite_BUNDLE=\"`echo \"$0\" | sed -e 's/\\/Contents\\/Resources\\/bin\\/startup.sh//'`\"\n")
    file(APPEND ${mac_startup} "NGSuite_RESOURCES=\"$NGSuite_BUNDLE/Contents/Resources\"\n")
    file(APPEND ${mac_startup} "export DYLD_LIBRARY_PATH=$NGSuite_RESOURCES/lib:$DYLD_LIBRARY_PATH\n")
    file(APPEND ${mac_startup} "export NETGENDIR=$NGSuite_RESOURCES/bin\n")
    file(APPEND ${mac_startup} "cd $NGSuite_RESOURCES\n") 
    file(APPEND ${mac_startup} "$NGSuite_RESOURCES/bin/netgen\n") 
    install(PROGRAMS ${mac_startup} DESTINATION bin)

    set(mac_ngsuite ${CMAKE_CURRENT_BINARY_DIR}/ngsuite.sh)
    file(WRITE ${mac_ngsuite} "#!/bin/sh\n")
    file(APPEND ${mac_ngsuite} "NGSuite_BUNDLE=\"`echo \"$0\" | sed -e 's/\\/Contents\\/MacOS\\/NGSuite//'`\"\n")
    file(APPEND ${mac_ngsuite} "NGSuite_RESOURCES=\"$NGSuite_BUNDLE/Contents/Resources\"\n")
    file(APPEND ${mac_ngsuite} "open -a /Applications/Utilities/Terminal.app $NGSuite_RESOURCES/bin/startup.sh\n")
    set(CPACK_BUNDLE_STARTUP_COMMAND ${mac_ngsuite})

    set(mac_plist ${CMAKE_CURRENT_BINARY_DIR}/Info.plist)
    file(WRITE ${mac_plist} "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n")
    file(APPEND ${mac_plist} "<!DOCTYPE plist PUBLIC \"-//Apple Computer//DTD PLIST 1.0//EN\" \"http://www.apple.com/DTDs/PropertyList-1.0.dtd\">\n")
    file(APPEND ${mac_plist} "<plist version=\"1.0\">\n")
    file(APPEND ${mac_plist} "<dict>\n")
    file(APPEND ${mac_plist} "	<key>CFBundleDevelopmentRegion</key>\n")
    file(APPEND ${mac_plist} "	<string>English</string>\n")
    file(APPEND ${mac_plist} "	<key>CFBundleExecutable</key>\n")
    file(APPEND ${mac_plist} "	<string>NGSuite</string>\n")
    file(APPEND ${mac_plist} "</dict>\n")
    file(APPEND ${mac_plist} "</plist>\n")
    set(CPACK_BUNDLE_PLIST ${mac_plist})
#    set(CPACK_DMG_FORMAT "UDRW")
#    set(CPACK_DMG_DS_STORE .DS_Store_dmg) 
#    set(CPACK_DMG_BACKGROUND_IMAGE sculpture312.png)

endif(APPLE)

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
string(TOUPPER "CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE}" name)
set(flags "${${name}} ${CMAKE_CXX_FLAGS}")
message("
------------------------------------------------------------------------
  ${PROJECT_NAME} ${PACKAGE_VERSION}:  Automatic configuration OK.

  Build type: ${CMAKE_BUILD_TYPE}
  Flags: ${flags}

  Enabled functionality:

    OCC: ............... ${USE_OCC}
    JPEGlib: ........... ${USE_JPEG}
    FFMPEG: ............ ${USE_MPEG}
    GUI: ............... ${USE_GUI}
    MPI: ............... ${USE_MPI}
    PYTHON: ............ ${USE_PYTHON}


  Building: ")
if(WIN32)
message("
    Open ${CMAKE_BINARY_DIR}/Netgen.sln and build solution to compile ${PROJECT_NAME}.

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

    ${CMAKE_INSTALL_PREFIX}

  Please set the following environment variables:

    NETGENDIR=${CMAKE_INSTALL_PREFIX}/bin")
if(USE_PYTHON)
message("    PYTHONPATH=.:${CMAKE_INSTALL_PREFIX}/${PYTHON_PACKAGES_INSTALL_DIR}")
endif(USE_PYTHON)
message("
------------------------------------------------------------------------
")

