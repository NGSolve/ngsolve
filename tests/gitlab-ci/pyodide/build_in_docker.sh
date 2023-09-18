#! /bin/bash

set -e

source /root/emsdk/emsdk_env.sh

export EMSCRIPTEN_SYSROOT=$(em-config CACHE)/sysroot
export EMSCRIPTEN_INCLUDE=$EMSCRIPTEN_SYSROOT/include
export EMSCRIPTEN_BIN=$EMSCRIPTEN_SYSROOT/bin
export EMSCRIPTEN_LIB=$EMSCRIPTEN_SYSROOT/lib/wasm32-emscripten/pic
export SIDE_MODULE_LDFLAGS="-O2 -g0 -sWASM_BIGINT -s SIDE_MODULE=1"
export TARGETINSTALLDIR=/root/xbuildenv/pyodide-root/cpython/installs/python-3.11.2/
export CMAKE_TOOLCHAIN_FILE=/root/emsdk/upstream/emscripten/cmake/Modules/Platform/Emscripten.cmake

sed -i 's/TARGET_SUPPORTS_SHARED_LIBS FALSE/TARGET_SUPPORTS_SHARED_LIBS TRUE/' $CMAKE_TOOLCHAIN_FILE
echo 'set(CMAKE_STRIP "${EMSCRIPTEN_ROOT_PATH}/emstrip${EMCC_SUFFIX}" CACHE FILEPATH "Emscripten strip")' >> ${CMAKE_TOOLCHAIN_FILE}
export CCACHE_DIR=/ccache
ccache -s

emcmake cmake /root/ngsolve \
      -DCMAKE_CROSSCOMPILING=ON \
      -DCMAKE_CXX_FLAGS="-sNO_DISABLE_EXCEPTION_CATCHING -I$TARGETINSTALLDIR/include/python3.11" \
      -DCMAKE_MODULE_LINKER_FLAGS="$SIDE_MODULE_LDFLAGS" \
      -DCMAKE_SHARED_LINKER_FLAGS="$SIDE_MODULE_LDFLAGS" \
      -DUSE_SUPERBUILD=ON \
      -DUSE_CCACHE=ON \
      -DNGLIB_LIBRARY_TYPE=STATIC \
      -DNGCORE_LIBRARY_TYPE=OBJECT \
      -DMAX_SYS_DIM=1 \
      -DCMAKE_BUILD_TYPE=Release \
      -DUSE_FASTCOMPILE=1 \
      -DUSE_CSG:UNINITIALIZED=ON \
      -DUSE_CUDA:BOOL=OFF \
      -DUSE_GUI:BOOL=OFF \
      -DUSE_HYPRE:BOOL=OFF \
      -DUSE_INTERFACE:UNINITIALIZED=ON \
      -DUSE_LAPACK:BOOL=ON \
      -DLAPACK_LIBRARIES="-lm" \
      -DUSE_MKL:BOOL=OFF \
      -DUSE_MPI:BOOL=OFF \
      -DUSE_MUMPS:BOOL=OFF \
      -DUSE_NATIVE_ARCH:UNINITIALIZED=OFF \
      -DUSE_PARDISO:BOOL=OFF \
      -DUSE_PYTHON:BOOL=ON \
      -DPython3_ROOT_DIR=$TARGETINSTALLDIR \
      -DPython3_LIBRARY=$TARGETINSTALLDIR/lib/libpython3.11.a \
      -DPython3_INCLUDE_DIR=$TARGETINSTALLDIR/include/python3.11 \
      -DPython3_INCLUDE_DIRS=$TARGETINSTALLDIR/include \
      -DNG_INSTALL_DIR_PYTHON=python \
      -DNG_INSTALL_DIR_LIB=python/netgen \
      -DUSE_STLGEOM:UNINITIALIZED=ON \
      -DUSE_UMFPACK:BOOL=OFF \
      -DNetgen_DIR=/opt/netgen/lib/cmake/netgen \
      -DCMAKE_INSTALL_PREFIX=/opt/netgen \
      -DCMAKE_POSITION_INDEPENDENT_CODE=ON \
      -DPYBIND_INCLUDE_DIR=/root/ngsolve/external_dependencies/netgen/external_dependencies/pybind11/include \
      -Dpybind11_INCLUDE_DIR=/root/ngsolve/external_dependencies/netgen/external_dependencies/pybind11/include \
      -DCMAKE_CXX_FLAGS_RELEASE="" \
      -DBUILD_FOR_CONDA=ON \
      -DUSE_OCC:BOOL=ON \
      -DOpenCascade_DIR=/opt/opencascade/lib/cmake/opencascade \
      -DBUILD_OCC=OFF \
      -DBUILD_ZLIB=ON \

emmake make -j9 install
