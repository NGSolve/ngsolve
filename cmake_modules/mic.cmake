set(CMAKE_SYSTEM_NAME Linux)
set(CMAKE_SYSTEM_PROCESSOR k1om)
set(CMAKE_SYSTEM_VERSION 1)

# specify the cross compiler
set(CMAKE_C_COMPILER   icc)
set(CMAKE_CXX_COMPILER icpc)
set(MPI_C_COMPILER mpiicc)
set(_CMAKE_TOOLCHAIN_PREFIX  x86_64-k1om-linux-)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mmic")
set(CMAKE_CXX_LINK_FLAGS "${CMAKE_CXX_LINK_FLAGS} -limf -lsvml -lirng -lintlc")

# where is the target environment 
set(CMAKE_FIND_ROOT_PATH /usr/linux-k1om-4.7)

