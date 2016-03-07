
if(APPLE)
    set(CMAKE_C_COMPILER "/usr/bin/gcc")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -x objective-c")
endif(APPLE)

if(WIN32)
  add_definitions("-DBUILD_togl -DUNICODE -D_UNICODE -DTOGL_USE_FONTS=0 -DSTDC_HEADERS -DSTDC_HEADER")
  add_library(togl SHARED togl.c toglProcAddr.c toglStubInit.c)
  install(TARGETS togl DESTINATION ${ng_install_dir} COMPONENT netgen)
else(WIN32)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fomit-frame-pointer -Wno-implicit-int")
  add_definitions("-DPACKAGE_NAME=\"Togl\" -DPACKAGE_TARNAME=\"togl\" -DPACKAGE_VERSION=\"2.1\" -DPACKAGE_STRING=\"Togl\ 2.1\" -DPACKAGE_BUGREPORT=\"\" -DPACKAGE_URL=\"\" -DSTDC_HEADERS=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=0 -DHAVE_LIMITS_H=1 -DHAVE_SYS_PARAM_H=1 -DUSE_THREAD_ALLOC=1 -D_REENTRANT=1 -D_THREAD_SAFE=1 -DTCL_THREADS=1 -DMODULE_SCOPE=extern\ __attribute__\(\(__visibility__\(\"hidden\"\)\)\) -D_LARGEFILE64_SOURCE=1 -DTCL_WIDE_INT_IS_LONG=1 -DUSE_TCL_STUBS=1 -DUSE_TK_STUBS=1 -DAUTOSTEREOD=\"\"") 

  include_directories(BEFORE "${TCL_INCLUDE_PATH}/tcl-private/generic" "${TCL_INCLUDE_PATH}/tcl-private/unix")
  include_directories(BEFORE "${TCL_INCLUDE_PATH}/tk-private/generic" "${TCL_INCLUDE_PATH}/tk-private/unix")
  include_directories(BEFORE "${TCL_INCLUDE_PATH}/tk-private/generic/ttk") 
  include_directories(BEFORE "${TK_INCLUDE_PATH}/../PrivateHeaders") 
  include_directories(BEFORE "${TCL_INCLUDE_PATH}") 
  include_directories(BEFORE "${TK_INCLUDE_PATH}") 

  add_library(togl togl.c toglProcAddr.c toglStubInit.c)
  target_link_libraries(togl -ldl)
endif(WIN32)

target_link_libraries(togl ${OPENGL_LIBRARIES} ${TK_LIBRARY} ${TCL_LIBRARY} )
set_target_properties(togl PROPERTIES POSITION_INDEPENDENT_CODE ON )
