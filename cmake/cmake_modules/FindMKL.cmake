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
#   MKL_LINK_FLAGS       : flags to use when linking with mkl
#   MKL_INCLUDE_DIRS     : set when MKL_INCLUDE_DIR found
#   MKL_LIBRARIES        : the library to link against.

# $(MKLROOT)/lib/intel64/libmkl_intel.a $(MKLROOT)/lib/intel64/libmkl_core.a $(MKLROOT)/lib/intel64/libmkl_gnu_thread.a -ldl -lpthread -lm
# $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_core.a $(MKLROOT)/lib/intel64/libmkl_gnu_thread.a -ldl -lpthread -lm
# $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_core.a $(MKLROOT)/lib/intel64/libmkl_gnu_thread.a -ldl -lpthread -lm

# interface layer, either lp64 or ilp64 (only for 64 bit build)
if(NOT MKL_INTERFACE)
    set(MKL_INTERFACE "lp64")
endif(NOT MKL_INTERFACE)

set(MKL_ARCH "intel64")

include(FindPackageHandleStandardArgs)

if(NOT MKL_LINK)
  set(MKL_LINK sdl)
endif()
set(d lib/cmake/mkl)
find_package(MKL CONFIG HINTS ${MKL_ROOT}/${d} PATHS /opt/intel/oneapi/mkl/latest/${d} ~/.local/${d})

if(NOT MKL_ROOT)
    set(MKL_ROOT $ENV{MKLROOT} CACHE PATH "Folder contains MKL")
endif(NOT MKL_ROOT)

if(MKL_FOUND)
  message("Found oneapi mkl")
  get_target_property(MKL_INCLUDE_DIRS MKL::MKL INTERFACE_INCLUDE_DIRECTORIES)
  set(MKL_LIBRARIES MKL::MKL)
  set(MKL_MINIMAL_LIBRARIES MKL::MKL)
else()
  # legacy search

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

  if(MKL_SDL AND NOT MKL_STATIC)
      find_library(MKL_LIBRARY mkl_rt
          PATHS ${MKL_ROOT}/lib/${MKL_ARCH}/)

      set(MKL_MINIMAL_LIBRARY ${MKL_LIBRARY})
      get_filename_component(MKL_LIB_DIR ${MKL_MINIMAL_LIBRARY} PATH)
      get_filename_component(MKL_LIB_NAME ${MKL_MINIMAL_LIBRARY} NAME_WE)
      string(SUBSTRING "${MKL_LIB_NAME}" 3 -1 MKL_LIB_NAME)
      set(MKL_LINK_FLAGS "-L${MKL_LIB_DIR} -l${MKL_LIB_NAME}") # CACHE string "MKL link flags")
  else()
      ######################### Interface layer #######################
      # win32 link advisor:
      # mkl_intel_lp64.lib mkl_core.lib mkl_intel_thread.lib libiomp5md.lib -ldl
      set(MKL_INTERFACE_LIBNAME "mkl_intel_${MKL_INTERFACE}")

      find_library(MKL_INTERFACE_LIBRARY ${MKL_INTERFACE_LIBNAME}
          PATHS ${MKL_ROOT}/lib/${MKL_ARCH}/)
      get_filename_component(MKL_LIBRARY_DIR ${MKL_INTERFACE_LIBRARY} PATH)

      ######################## Threading layer ########################
      if(MKL_MULTI_THREADED)
        if(WIN32)
          set(MKL_THREADING_LIBNAME mkl_intel_thread)
        else(WIN32)
          set(MKL_THREADING_LIBNAME mkl_gnu_thread)
        endif()
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
      find_library(MKL_SCALAPACK_LIBRARY mkl_scalapack_${MKL_INTERFACE}
          PATHS ${MKL_ROOT}/lib/${MKL_ARCH}/)

      find_library(MKL_BLACS_LIBRARY mkl_blacs_openmpi_${MKL_INTERFACE}
          PATHS ${MKL_ROOT}/lib/${MKL_ARCH}/)

  #     ############################ RTL layer ##########################
      if(WIN32)
          set(MKL_RTL_LIBNAME libiomp5md)
          find_library(MKL_RTL_LIBRARY ${MKL_RTL_LIBNAME} PATHS ${MKL_ROOT}/../compiler/lib/${MKL_ARCH})
          find_file(MKL_RTL_LIBRARY_BIN ${MKL_RTL_LIBNAME}.dll PATHS ${MKL_ROOT}/../redist/${MKL_ARCH}/compiler)
      endif()
  #         set(MKL_RTL_LIBNAME libiomp5)
  #     endif()
  #     find_library(MKL_RTL_LIBRARY ${MKL_RTL_LIBNAME}
  #         PATHS ${INTEL_RTL_ROOT}/lib)

      if(WIN32)
          set(MKL_MINIMAL_LIBRARY ${MKL_INTERFACE_LIBRARY} ${MKL_THREADING_LIBRARY} ${MKL_CORE_LIBRARY} ${MKL_RTL_LIBRARY})
          set(MKL_LIBRARY ${MKL_INTERFACE_LIBRARY} ${MKL_THREADING_LIBRARY} ${MKL_CORE_LIBRARY} ${MKL_RTL_LIBRARY})
      else(WIN32)
          set(MKL_LIBRARY "-Wl,--start-group ${MKL_INTERFACE_LIBRARY} ${MKL_THREADING_LIBRARY} ${MKL_CORE_LIBRARY} ${MKL_FFT_LIBRARY} ${MKL_BLACS_LIBRARY} ${MKL_SCALAPACK_LIBRARY} ${MKL_RTL_LIBRARY} -Wl,--end-group -ldl -lpthread -lm -fopenmp")
          set(MKL_LINK_FLAGS "-Wl,--start-group ${MKL_INTERFACE_LIBRARY} ${MKL_THREADING_LIBRARY} ${MKL_CORE_LIBRARY} -Wl,--end-group -ldl -lpthread -lm -fopenmp" )
          set(MKL_MINIMAL_LIBRARY ${MKL_LINK_FLAGS})
  #  ${MKL_ROOT}/lib/${MKL_ARCH}/libmkl_blacs_openmpi_lp64.a 
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
endif()
