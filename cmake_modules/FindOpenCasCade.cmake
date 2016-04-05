# Try to find OCC
# Once done this will define
#
# OCC_FOUND          - system has OCC - OpenCASCADE
# OCC_INCLUDE_DIR    - where the OCC include directory can be found
# OCC_LIBRARY_DIR    - where the OCC library directory can be found
# OCC_LIBRARIES      - Link this to use OCC

if(WIN32)
    find_path(OCC_INCLUDE_DIR Standard_Version.hxx PATH_SUFFIXES inc ../inc)
    find_library(OCC_LIBRARY TKernel)
else(WIN32)
    find_path(OCC_INCLUDE_DIR Standard_Version.hxx
      /usr/include/opencascade
      /usr/local/include/opencascade
      /usr/include/oce
      /usr/local/include/oce
      /opt/opencascade/include
      /opt/opencascade/inc
    )
    find_library(OCC_LIBRARY TKernel
      /usr/lib
      /usr/local/lib
      /opt/opencascade/lib
    )
endif(WIN32)

if(OCC_LIBRARY)
    get_filename_component(OCC_LIBRARY_DIR ${OCC_LIBRARY} PATH)
endif(OCC_LIBRARY)

if(OCC_INCLUDE_DIR)
    file(STRINGS ${OCC_INCLUDE_DIR}/Standard_Version.hxx OCC_MAJOR
      REGEX "#define OCC_VERSION_MAJOR.*"
    )
    string(REGEX MATCH "[0-9]+" OCC_MAJOR ${OCC_MAJOR})
    file(STRINGS ${OCC_INCLUDE_DIR}/Standard_Version.hxx OCC_MINOR
      REGEX "#define OCC_VERSION_MINOR.*"
    )
    string(REGEX MATCH "[0-9]+" OCC_MINOR ${OCC_MINOR})
    file(STRINGS ${OCC_INCLUDE_DIR}/Standard_Version.hxx OCC_MAINT
      REGEX "#define OCC_VERSION_MAINTENANCE.*"
    )
    string(REGEX MATCH "[0-9]+" OCC_MAINT ${OCC_MAINT})

    set(OCC_VERSION_STRING "${OCC_MAJOR}.${OCC_MINOR}.${OCC_MAINT}")
endif(OCC_INCLUDE_DIR)


set(OCC_LIBRARY_NAMES
    TKBO
    TKBool
    TKBRep
    TKCAF
    TKCDF
    TKernel
    TKG2d
    TKG3d
    TKGeomAlgo
    TKGeomBase
    TKHLR
    TKIGES
    TKLCAF
    TKMath
    TKMesh
    TKOffset
    TKPrim
    TKService
    TKShHealing
    TKSTEP
    TKSTEP209
    TKSTEPAttr
    TKSTEPBase
    TKSTL
    TKTopAlgo
    TKV3d
    TKXCAF
    TKXDEIGES
    TKXDESTEP
    TKXSBase
)

foreach( libname ${OCC_LIBRARY_NAMES} )
    find_library( ${libname} ${libname} ${OCC_LIBRARY_DIR} )
    set(OCC_LIBRARIES ${OCC_LIBRARIES} ${${libname}})
endforeach()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(OCC REQUIRED_VARS OCC_INCLUDE_DIR VERSION_VAR OCC_VERSION_STRING ${OCC_LIBRARIY_NAMES})

if(OCC_FOUND)
    message(STATUS "-- Found OpenCASCADE version: ${OCC_VERSION_STRING}")
    message(STATUS "-- OpenCASCADE include directory: ${OCC_INCLUDE_DIR}")
    message(STATUS "-- OpenCASCADE shared libraries directory: ${OCC_LIBRARY_DIR}")
    message(STATUS "-- OpenCASCADE shared libraries :\n ${OCC_LIBRARIES}")
endif(OCC_FOUND)

