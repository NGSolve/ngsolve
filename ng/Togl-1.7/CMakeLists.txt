add_definitions("-DPACKAGE_NAME=\"Togl\" -DPACKAGE_TARNAME=\"togl\" -DPACKAGE_VERSION=\"1.7\" -DPACKAGE_STRING=\"Togl\ 1.7\" -DPACKAGE_BUGREPORT=\"\" -DSTDC_HEADERS=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_LIMITS_H=1 -DHAVE_SYS_PARAM_H=1 -DUSE_THREAD_ALLOC=1 -D_REENTRANT=1 -D_THREAD_SAFE=1 -DTCL_THREADS=1 -D_LARGEFILE64_SOURCE=1 -DTCL_WIDE_INT_IS_LONG=1 -DUSE_TCL_STUBS=1 -DUSE_TK_STUBS=1")

# include_directories("/usr/include/tcl8.5" "/usr/include/tcl8.5/tk-private/generic" "/usr/include/tcl8.5/tk-private/unix")
# SET(CMAKE_CXX_FLAGS "-O2 -fomit-frame-pointer -Wall -Wno-implicit-int -fPIC -c")
include_directories("${TCL_INCLUDE_PATH}/tk-private/generic" "${TCL_INCLUDE_PATH}/tk-private/unix")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -fomit-frame-pointer -Wall -Wno-implicit-int")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O2 -fomit-frame-pointer -Wall -Wno-implicit-int")
add_library(togl togl.c)
target_link_libraries(togl ${OPENGL_LIBRARIES})
set_target_properties(togl PROPERTIES POSITION_INDEPENDENT_CODE ON )
# 
# gcc -DPACKAGE_NAME=\"Togl\" -DPACKAGE_TARNAME=\"togl\" -DPACKAGE_VERSION=\"1.7\" -DPACKAGE_STRING=\"Togl\ 1.7\" -DPACKAGE_BUGREPORT=\"\" -DSTDC_HEADERS=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_LIMITS_H=1 -DHAVE_SYS_PARAM_H=1 -DUSE_THREAD_ALLOC=1 -D_REENTRANT=1 -D_THREAD_SAFE=1 -DTCL_THREADS=1 -D_LARGEFILE64_SOURCE=1 -DTCL_WIDE_INT_IS_LONG=1 -DUSE_TCL_STUBS=1 -DUSE_TK_STUBS=1                                                                                                                         
# -I"/usr/include/tcl8.6" -I"/usr/include/tcl8.6/tk-private/generic" -I"/usr/include/tcl8.6/tk-private/unix" 
# -O2 -fomit-frame-pointer -Wall -Wno-implicit-int -fPIC -c `echo togl.c` -o togl.o
# rm -f libTogl1.7.so
# gcc -pipe -shared -o libTogl1.7.so togl.o -lX11 -lGL -lXmu  -L/usr/lib/x86_64-linux-gnu -ltclstub8.6 -L/usr/lib/x86_64-linux-gnu -ltkstub8.6 
# : libTogl1.7.so


