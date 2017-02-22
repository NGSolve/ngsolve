set(CPACK_INSTALL_CMAKE_PROJECTS 
  "${CMAKE_BINARY_DIR}/netgen/netgen;Netgen;ALL;/"
  "${CMAKE_BINARY_DIR}/ngsolve;NGSolve;ALL;/"
)

if(NGSOLVE_VERSION_TWEAK)
  # nightly build
  set(CPACK_PACKAGE_VERSION ${NGSOLVE_VERSION_LONG})
  set(CPACK_SOURCE_PACKAGE_FILE_NAME ngsolve_v${NGSOLVE_VERSION_LONG})
  set(CPACK_DEBIAN_PACKAGE_NAME ngsolve_nightly)
else(NGSOLVE_VERSION_TWEAK)
  # release
  set(CPACK_PACKAGE_NAME ngsolve)
  set(CPACK_DEBIAN_PACKAGE_NAME ngsolve)
  set(CPACK_PACKAGE_VERSION ${NGSOLVE_VERSION})
  set(CPACK_SOURCE_PACKAGE_FILE_NAME ngsolve_v${NGSOLVE_VERSION})
endif(NGSOLVE_VERSION_TWEAK)


# set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "automatic 3d tetrahedral mesh generator")
# set(CPACK_PACKAGE_DESCRIPTION "NETGEN is an automatic 3d tetrahedral mesh generator. It accepts input from constructive solid geometry (CSG) or boundary representation (BRep) from STL file format. The connection to a geometry kernel allows the handling of IGES and STEP files. NETGEN contains modules for mesh optimization and hierarchical mesh refinement. Netgen is open source based on the LGPL license. It is available for Unix/Linux and Windows.")

set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "NGSolve Finite Element Library")
set(CPACK_PACKAGE_DESCRIPTION "NGSolve is a general purpose Finite Element Library on top of Netgen. With the basic library one can solve heat flow equations, Maxwell equations, and solid mechanical problems. Several add-ons are available for particular application classes.") 


set(CPACK_COMPONENTS_ALL netgen netgen_devel netgen_tutorial ngsolve ngsolve_devel ngsolve_tutorial)
set(CPACK_PACKAGE_VENDOR "Vienna University of Technology")

set(CPACK_SOURCE_GENERATOR "TGZ")
set(CPACK_SOURCE_IGNORE_FILES "/.dmg/;/.msi/;/.deb/;/.git/;/cmake/;/build/;/.gz/;/.zip/;~$;${CPACK_SOURCE_IGNORE_FILES}")

if(WIN32)
  include(cmake_modules/package_windows.cmake)
elseif(APPLE)
  include(cmake_modules/package_macos.cmake)
elseif(UNIX)
  include(cmake_modules/package_linux.cmake)
endif()

