set(CPACK_INSTALL_CMAKE_PROJECTS 
  "${CMAKE_BINARY_DIR}/netgen/netgen;Netgen;ALL;/"
  "${CMAKE_BINARY_DIR}/ngsolve;NGSolve;ALL;/"
)

set(CPACK_PACKAGE_NAME NGSuite)
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "automatic 3d tetrahedral mesh generator")
set(CPACK_PACKAGE_DESCRIPTION "NETGEN is an automatic 3d tetrahedral mesh generator. It accepts input from constructive solid geometry (CSG) or boundary representation (BRep) from STL file format. The connection to a geometry kernel allows the handling of IGES and STEP files. NETGEN contains modules for mesh optimization and hierarchical mesh refinement. Netgen is open source based on the LGPL license. It is available for Unix/Linux and Windows.")

set(CPACK_COMPONENTS_ALL netgen netgen_devel netgen_tutorial ngsolve ngsolve_devel ngsolve_tutorial)
set(CPACK_PACKAGE_VENDOR "Vienna University of Technology")

set(CPACK_SOURCE_GENERATOR "TGZ")
set(CPACK_SOURCE_IGNORE_FILES "/cmake/;/build/;/.gz/;/.zip/;~$;${CPACK_SOURCE_IGNORE_FILES}")

if(WIN32)
  include(cmake_modules/package_windows.cmake)
elseif(APPLE)
  include(cmake_modules/package_macos.cmake)
elseif(UNIX)
  include(cmake_modules/package_linux.cmake)
endif()

include(CPack)
