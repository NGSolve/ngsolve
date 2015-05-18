# - Find Intel MKL
# Find the MKL libraries
#
# Options:
#
#   MKL_STATIC       :   use static linking
#   MKL_MULTI_THREADED:   use multi-threading
#   MKL_SDL           :   Single Dynamic Library interface
#   MKL_ARCH : architecture to link (possible values: intel64, mic)
#
# This module defines the following variables:
#
#   MKL_FOUND            : True if MKL_INCLUDE_DIR are found
#   MKL_INCLUDE_DIR      : where to find mkl.h, etc.
#   MKL_INCLUDE_DIRS     : set when MKL_INCLUDE_DIR found
#   MKL_LIBRARIES        : the library to link against.

# $(MKLROOT)/lib/intel64/libmkl_intel.a $(MKLROOT)/lib/intel64/libmkl_core.a $(MKLROOT)/lib/intel64/libmkl_gnu_thread.a -ldl -lpthread -lm
# $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_core.a $(MKLROOT)/lib/intel64/libmkl_gnu_thread.a -ldl -lpthread -lm
# $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_core.a $(MKLROOT)/lib/intel64/libmkl_gnu_thread.a -ldl -lpthread -lm

# interface layer, either _lp64 or _lp64
set(MKL_INTERFACE_LAYER "_lp64")

if(NOT MKL_ARCH)
    set(MKL_ARCH "intel64")
endif(NOT MKL_ARCH)

include(FindPackageHandleStandardArgs)

if(NOT MKL_ROOT)
    set(MKL_ROOT $ENV{MKLROOT} CACHE PATH "Folder contains MKL")
endif(NOT MKL_ROOT)

# Find include dir
find_path(MKL_INCLUDE_DIR mkl.h PATHS ${MKL_ROOT}/include)

# Find libraries

# Handle suffix
set(_MKL_ORIG_CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES})

if(WIN32)
    if(MKL_STATIC)
        set(CMAKE_FIND_LIBRARY_SUFFIXES .lib)
    elseif(MKL_SDL)
        set(CMAKE_FIND_LIBRARY_SUFFIXES .lib)
    else()
        set(CMAKE_FIND_LIBRARY_SUFFIXES _dll.lib)
    endif()
elseif(APPLE)
    if(MKL_STATIC)
        set(CMAKE_FIND_LIBRARY_SUFFIXES .a)
    else()
        set(CMAKE_FIND_LIBRARY_SUFFIXES .dylib)
    endif()
else()
    if(MKL_STATIC)
        set(CMAKE_FIND_LIBRARY_SUFFIXES .a)
    else()
        set(CMAKE_FIND_LIBRARY_SUFFIXES .so)
    endif()
endif()


# MKL is composed by four layers: Interface, Threading, Computational and RTL

if(MKL_SDL)
    find_library(MKL_LIBRARY mkl_rt
        PATHS ${MKL_ROOT}/lib/${MKL_ARCH}/)

    set(MKL_MINIMAL_LIBRARY ${MKL_LIBRARY})
else()
    ######################### Interface layer #######################
    # win32 link advisor:
    # mkl_intel_lp64.lib mkl_core.lib mkl_intel_thread.lib libiomp5md.lib -ldl
    set(MKL_INTERFACE_LIBNAME "mkl_intel${MKL_INTERFACE_LAYER}")

    find_library(MKL_INTERFACE_LIBRARY ${MKL_INTERFACE_LIBNAME}
        PATHS ${MKL_ROOT}/lib/${MKL_ARCH}/)

    ######################## Threading layer ########################
    if(MKL_MULTI_THREADED)
        if(WIN32)
            set(MKL_THREADING_LIBNAME mkl_intel_thread)
        else(WIN32)
            set(MKL_THREADING_LIBNAME mkl_gnu_thread)
        endif(WIN32)
    else()
        set(MKL_THREADING_LIBNAME mkl_sequential)
    endif()

    find_library(MKL_THREADING_LIBRARY ${MKL_THREADING_LIBNAME}
        PATHS ${MKL_ROOT}/lib/${MKL_ARCH}/)

    ####################### Computational layer #####################
    find_library(MKL_CORE_LIBRARY mkl_core
        PATHS ${MKL_ROOT}/lib/${MKL_ARCH}/)
    find_library(MKL_FFT_LIBRARY mkl_cdft_core
        PATHS ${MKL_ROOT}/lib/${MKL_ARCH}/)
    find_library(MKL_SCALAPACK_LIBRARY mkl_scalapack${MKL_INTERFACE_LAYER}
        PATHS ${MKL_ROOT}/lib/${MKL_ARCH}/)

#     ############################ RTL layer ##########################
    if(WIN32)
        set(MKL_RTL_LIBNAME libiomp5md)
        find_library(MKL_RTL_LIBRARY ${MKL_RTL_LIBNAME}
            PATHS ${MKL_ROOT}/../compiler/lib/${MKL_ARCH})
    endif()
#         set(MKL_RTL_LIBNAME libiomp5)
#     endif()
#     find_library(MKL_RTL_LIBRARY ${MKL_RTL_LIBNAME}
#         PATHS ${INTEL_RTL_ROOT}/lib)

    if(WIN32)
        set(MKL_MINIMAL_LIBRARY ${MKL_INTERFACE_LIBRARY} ${MKL_THREADING_LIBRARY} ${MKL_CORE_LIBRARY} ${MKL_RTL_LIBRARY})
        set(MKL_LIBRARY ${MKL_INTERFACE_LIBRARY} ${MKL_THREADING_LIBRARY} ${MKL_CORE_LIBRARY} ${MKL_RTL_LIBRARY})
    else(WIN32)
        set(MKL_LIBRARY "-Wl,--start-group ${MKL_INTERFACE_LIBRARY} ${MKL_THREADING_LIBRARY} ${MKL_CORE_LIBRARY} ${MKL_FFT_LIBRARY} ${MKL_SCALAPACK_LIBRARY} ${MKL_RTL_LIBRARY} -Wl,--end-group -ldl -lpthread -lm")
        set(MKL_MINIMAL_LIBRARY "-Wl,--start-group ${MKL_INTERFACE_LIBRARY} ${MKL_THREADING_LIBRARY} ${MKL_CORE_LIBRARY} -Wl,--end-group ${MKL_ROOT}/lib/${MKL_ARCH}/libmkl_blacs_openmpi_lp64.a -ldl -lpthread -lm -fopenmp")
    endif(WIN32)
endif()

set(CMAKE_FIND_LIBRARY_SUFFIXES ${_MKL_ORIG_CMAKE_FIND_LIBRARY_SUFFIXES})

find_package_handle_standard_args(MKL DEFAULT_MSG
    MKL_INCLUDE_DIR MKL_LIBRARY MKL_MINIMAL_LIBRARY)

if(MKL_FOUND)
    set(MKL_INCLUDE_DIRS ${MKL_INCLUDE_DIR})
    set(MKL_LIBRARIES ${MKL_LIBRARY})
    set(MKL_MINIMAL_LIBRARIES ${MKL_MINIMAL_LIBRARY})
#     message("***************************************************MKL LIBS:\n${MKL_LIBRARIES}")
endif()
