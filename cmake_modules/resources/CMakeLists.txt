cmake_minimum_required(VERSION 2.8)
project(NGSuite)
set(NGSUITE_VERSION_MAJOR 6)
set(NGSUITE_VERSION_MINOR 1)
set(NGSUITE_VERSION "${NGSUITE_VERSION_MAJOR}.${NGSUITE_VERSION_MINOR}")
string(TIMESTAMP NGSUITE_VERSION_PATCH "%y%U%w" )
# set(NGSUITE_VERSION_PATCH 15425)
# set(NGSUITE_VERSION_LONG "${NGSUITE_VERSION_MAJOR}.${NGSUITE_VERSION_MINOR}.${NGSUITE_VERSION_PATCH}")
set(NGSUITE_VERSION_LONG "${NGSUITE_VERSION}-dev")
set(PACKAGE_VERSION "${NGSUITE_VERSION}.${NGSUITE_VERSION_PATCH}")
include(ExternalProject)

option( USE_OCC     "compile with OpenCascade geometry kernel" OFF)

if(NOT NETGEN_GIT_BRANCH)
    set(NETGEN_GIT_BRANCH master)
endif(NOT NETGEN_GIT_BRANCH)

if(NOT NGSOLVE_GIT_BRANCH)
    set(NGSOLVE_GIT_BRANCH master)
endif(NOT NGSOLVE_GIT_BRANCH)

# set(INSTALL_DIR ${CMAKE_BINARY_DIR}/install)

list(APPEND COMMON_CMAKE_ARGS "-DINSTALL_DIR=${CMAKE_BINARY_DIR}/install" "-DUSE_NATIVE_ARCH=OFF" "-DCMAKE_CXX_FLAGS_RELEASE=\"/W0\"")
list(APPEND NETGEN_CMAKE_ARGS "-DUSE_OCC=${USE_OCC}")
list(APPEND NGSOLVE_CMAKE_ARGS "-DNETGEN_SOURCE_DIR=${CMAKE_BINARY_DIR}/Netgen-prefix/src/Netgen" "-DUSE_MKL=ON" "-DMKL_STATIC=ON" "-DUSE_CPP_THREADS=ON")

foreach( arg IN LISTS COMMON_CMAKE_ARGS NETGEN_CMAKE_ARGS )
  list(APPEND NETGEN_ARGS "${arg}")
endforeach()

foreach( arg IN LISTS COMMON_CMAKE_ARGS NGSOLVE_CMAKE_ARGS )
  list(APPEND NGSOLVE_ARGS "${arg}")
endforeach()

if(WIN32)
else(WIN32)
    set(NGSOLVE_CMAKE_ARGS ${NGSOLVE_CMAKE_ARGS} -DMKL_ROOT=/opt/intel/mkl) 
endif(WIN32)

file(TO_CMAKE_PATH "${ADDITIONAL_PATHS}" ADDITIONAL_PATHS)

message("ARGS\n${NETGEN_ARGS}")
ExternalProject_Add(
    Netgen

#     SOURCE_DIR "C:/Users/win7_work/src/netgen/"
#     DOWNLOAD_COMMAND ""
#     UPDATE_COMMAND ""

#     GIT_REPOSITORY git://git.code.sf.net/p/netgen-mesher/git

    UPDATE_COMMAND "git" "pull"
    GIT_REPOSITORY git@data.asc.tuwien.ac.at:jschoeberl/netgen.git
    GIT_TAG ${NETGEN_GIT_BRANCH}
    CMAKE_GENERATOR ${CMAKE_GENERATOR}
    CMAKE_ARGS ${NETGEN_ARGS}
    CMAKE_CACHE_ARGS "-DADDITIONAL_PATHS:PATH=${ADDITIONAL_PATHS}"
)

ExternalProject_Add(
    NGSolve
    DEPENDS Netgen

#     SOURCE_DIR "C:/Users/Matthias Hochsteger/src/ngsolve/"
#     DOWNLOAD_COMMAND ""
#     UPDATE_COMMAND ""

#     GIT_REPOSITORY git://git.code.sf.net/p/ngsolve/git

    UPDATE_COMMAND "git" "pull"
    GIT_REPOSITORY git@data.asc.tuwien.ac.at:jschoeberl/ngsolve.git
    GIT_TAG ${NGSOLVE_GIT_BRANCH}
    CMAKE_GENERATOR ${CMAKE_GENERATOR}
    CMAKE_GENERATOR_TOOLSET ${CMAKE_GENERATOR_TOOLSET}
    CMAKE_ARGS ${NGSOLVE_ARGS}
    CMAKE_CACHE_ARGS "-DADDITIONAL_PATHS:PATH=${ADDITIONAL_PATHS}"
) 

set(CPACK_INSTALL_CMAKE_PROJECTS 
    "${CMAKE_BINARY_DIR}/Netgen-prefix/src/Netgen-build;Netgen;ALL;/"
    "${CMAKE_BINARY_DIR}/NGSolve-prefix/src/NGSolve-build;NGSolve;ALL;/"
)


#######################################################################
# CPACK configuration

set(CPACK_COMPONENTS_ALL netgen netgen_devel netgen_tutorial ngsolve ngsolve_devel ngsolve_tutorial)

set(CPACK_COMPONENT_GROUP_GROUP_NETGEN_DISPLAY_NAME "Netgen")
set(CPACK_COMPONENT_NETGEN_GROUP "group_netgen")
set(CPACK_COMPONENT_NETGEN_DEVEL_GROUP "group_netgen")
set(CPACK_COMPONENT_NETGEN_TUTORIAL_GROUP "group_netgen")

set(CPACK_COMPONENT_NETGEN_DISPLAY_NAME "Runtime")
set(CPACK_COMPONENT_NETGEN_DEVEL_DISPLAY_NAME "Development files")
set(CPACK_COMPONENT_NETGEN_TUTORIAL_DISPLAY_NAME "Tutorials")

set(CPACK_COMPONENT_NETGEN_DEVEL_DEPENDS netgen)
set(CPACK_COMPONENT_NETGEN_TUTORIAL_DEPENDS netgen)

set(CPACK_COMPONENT_GROUP_GROUP_NGSOLVE_DISPLAY_NAME "NGSolve")
set(CPACK_COMPONENT_NGSOLVE_GROUP "group_ngsolve")
set(CPACK_COMPONENT_NGSOLVE_DEVEL_GROUP "group_ngsolve")
set(CPACK_COMPONENT_NGSOLVE_TUTORIAL_GROUP "group_ngsolve")

set(CPACK_COMPONENT_NGSOLVE_DISPLAY_NAME "Runtime")
set(CPACK_COMPONENT_NGSOLVE_DEVEL_DISPLAY_NAME "Development files")
set(CPACK_COMPONENT_NGSOLVE_TUTORIAL_DISPLAY_NAME "Tutorials")

set(CPACK_COMPONENT_NGSOLVE_DEPENDS netgen)
set(CPACK_COMPONENT_NGSOLVE_DEVEL_DEPENDS ngsolve)
set(CPACK_COMPONENT_NGSOLVE_TUTORIAL_DEPENDS ngsolve)


set(CPACK_SOURCE_GENERATOR "TGZ")
set(CPACK_SOURCE_IGNORE_FILES "/cmake/;/build/;/.gz/;~$;${CPACK_SOURCE_IGNORE_FILES}")
# set(CPACK_PACKAGE_VERSION ${PACKAGE_VERSION})
set(CPACK_PACKAGE_VERSION_MAJOR ${NGSUITE_VERSION_MAJOR})
set(CPACK_PACKAGE_VERSION_MINOR ${NGSUITE_VERSION_MINOR})
set(CPACK_PACKAGE_VERSION_PATCH ${NGSUITE_VERSION_PATCH})
# set(CPACK_PACKAGING_INSTALL_PREFIX "${INSTALL_DIR}")

set(CPACK_PACKAGE_NAME NGSuite)
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "automatic 3d tetrahedral mesh generator")
set(CPACK_PACKAGE_DESCRIPTION "NETGEN is an automatic 3d tetrahedral mesh generator. It accepts input from constructive solid geometry (CSG) or boundary representation (BRep) from STL file format. The connection to a geometry kernel allows the handling of IGES and STEP files. NETGEN contains modules for mesh optimization and hierarchical mesh refinement. Netgen is open source based on the LGPL license. It is available for Unix/Linux and Windows.")
if(UNIX)
    # Debian packager
    execute_process(COMMAND grep CODENAME /etc/lsb-release OUTPUT_VARIABLE temp OUTPUT_STRIP_TRAILING_WHITESPACE)
    if(temp)
        set(CPACK_GENERATOR "DEB")
        string(SUBSTRING ${temp} 17 -1 UBUNTU_VERSION)
        message("ubuntu version: ${UBUNTU_VERSION}")
        if(${UBUNTU_VERSION} STREQUAL "utopic")
            set(CPACK_DEBIAN_PACKAGE_DEPENDS "python3, libboost-python1.55.0, libtk8.5, libtcl8.5, tix, libxmu6")
        endif(${UBUNTU_VERSION} STREQUAL "utopic")
        if(${UBUNTU_VERSION} STREQUAL "trusty")
            set(CPACK_DEBIAN_PACKAGE_DEPENDS "python3, libboost-python1.54.0, libtk8.5, libtcl8.5, tix, libxmu6")
        endif(${UBUNTU_VERSION} STREQUAL "trusty")
        execute_process(COMMAND dpkg --print-architecture OUTPUT_VARIABLE CPACK_DEBIAN_PACKAGE_ARCHITECTURE OUTPUT_STRIP_TRAILING_WHITESPACE)
        set(CPACK_DEBIAN_PACKAGE_MAINTAINER Matthias Hochsteger <matthias.hochsteger@tuwien.ac.at>)
# TODO (use dpkg -s <packagename> to retrieve version)
        if(USE_MPI)
            set(CPACK_DEBIAN_PACKAGE_DEPENDS "${CPACK_DEBIAN_PACKAGE_DEPENDS}, libmetis5, openmpi-bin")
            set(CPACK_PACKAGE_NAME "${CPACK_PACKAGE_NAME}_mpi")
        endif(USE_MPI)
        set(CPACK_DEBIAN_PACKAGE_SECTION Science)
        set(CPACK_DEBIAN_PACKAGE_VERSION ${CPACK_PACKAGE_VERSION})
    endif(temp)
endif(UNIX)
if(WIN32)
    # WIX packager
    set(CPACK_WIX_PRODUCT_GUID 05BC95C0-440F-4D34-81DF-F54F966EC51E)
    set(CPACK_WIX_UPGRADE_GUID 5618ABB2-B105-4FEA-9A01-0E1094DC91FF)

    set(CPACK_GENERATOR "WIX")
    set(CPACK_WIX_PROGRAM_MENU_FOLDER "Netgen${NGSUITE_VERSION}-dev" )
    set(CPACK_PACKAGE_INSTALL_DIRECTORY "netgen${NGSUITE_VERSION}")
    set(CPACK_PACKAGE_VENDOR "Vienna University of Technology")
    set(CPACK_WIX_PATCH_FILE ${CMAKE_SOURCE_DIR}/wix_patch.xml)
    set(CPACK_RESOURCE_FILE_LICENSE ${CMAKE_SOURCE_DIR}/LICENSE.txt)
    set(CPACK_PACKAGE_EXECUTABLES ${CPACK_PACKAGE_EXECUTABLES} "netgen" "Netgen")
    set(CPACK_WIX_PRODUCT_ICON ${CMAKE_SOURCE_DIR}/netgen.ico)
    set(CPACK_WIX_UI_BANNER ${CMAKE_SOURCE_DIR}/banner.bmp)
    set(CPACK_WIX_UI_DIALOG ${CMAKE_SOURCE_DIR}/sculpture312_50.bmp)
    set(CPACK_WIX_TEMPLATE ${CMAKE_SOURCE_DIR}/WIX.template.in)
    set(CPACK_PACKAGE_FILE_NAME "NGSuite-${NGSUITE_VERSION}")
#     if(USE_OCC)
#       set(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_FILE_NAME}-OCC")
#     endif(USE_OCC)
    if( CMAKE_SIZEOF_VOID_P EQUAL 8)
        set(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_FILE_NAME}-win64")
    else( CMAKE_SIZEOF_VOID_P EQUAL 8)
        set(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_FILE_NAME}-win32")
    endif( CMAKE_SIZEOF_VOID_P EQUAL 8)

#     set(CPACK_WIX_EXTRA_SOURCES ${CMAKE_SOURCE_DIR}/extra_features.wxs)

endif(WIN32)

include(CPack)
message("CPACK_TOPLEVEL_TAG:  ${CPACK_TOPLEVEL_TAG}")
message("CPACK_SYSTEM_NAME:  ${CPACK_SYSTEM_NAME}")
message("CMAKE_SYSTEM_NAME:  ${CMAKE_SYSTEM_NAME}")

message("CMAKE_HOST_APPLE: ${CMAKE_HOST_APPLE}")
message("CMAKE_HOST_SYSTEM_NAME: ${CMAKE_HOST_SYSTEM_NAME}")
message("CMAKE_HOST_SYSTEM_PROCESSOR: ${CMAKE_HOST_SYSTEM_PROCESSOR}")
message("CMAKE_HOST_SYSTEM: ${CMAKE_HOST_SYSTEM}")
message("CMAKE_HOST_SYSTEM_VERSION: ${CMAKE_HOST_SYSTEM_VERSION}")
message("CMAKE_HOST_UNIX: ${CMAKE_HOST_UNIX}")
message("CMAKE_HOST_WIN32: ${CMAKE_HOST_WIN32}")
message("CMAKE_LIBRARY_ARCHITECTURE_REGEX: ${CMAKE_LIBRARY_ARCHITECTURE_REGEX}")
message("CMAKE_LIBRARY_ARCHITECTURE: ${CMAKE_LIBRARY_ARCHITECTURE}")
message("CMAKE_OBJECT_PATH_MAX: ${CMAKE_OBJECT_PATH_MAX}")
message("CMAKE_SYSTEM_NAME: ${CMAKE_SYSTEM_NAME}")
message("CMAKE_SYSTEM_PROCESSOR: ${CMAKE_SYSTEM_PROCESSOR}")
message("CMAKE_SYSTEM: ${CMAKE_SYSTEM}")
message("CMAKE_SYSTEM_VERSION: ${CMAKE_SYSTEM_VERSION}")


