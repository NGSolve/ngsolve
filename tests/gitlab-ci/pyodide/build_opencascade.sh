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

cd
wget https://github.com/Open-Cascade-SAS/OCCT/archive/refs/tags/V7_6_3.zip
unzip V7_6_3.zip
cd OCCT-7_6_3
mkdir build
cd build

emcmake cmake .. \
      -DCMAKE_CROSSCOMPILING=ON \
      -DCMAKE_CXX_FLAGS="-sNO_DISABLE_EXCEPTION_CATCHING -I$TARGETINSTALLDIR/include/python3.11" \
      -DCMAKE_INSTALL_PREFIX=/opt/opencascade \
      -DBUILD_LIBRARY_TYPE:STRING=Static \
      -DBUILD_MODULE_FoundationClasses:BOOL=ON \
      -DBUILD_MODULE_ModelingData:BOOL=ON \
      -DBUILD_MODULE_ModelingAlgorithms:BOOL=ON \
      -DBUILD_MODULE_DataExchange:BOOL=ON \
      -DBUILD_MODULE_Visualization:BOOL=OFF \
      -DBUILD_MODULE_ApplicationFramework:BOOL=OFF \
      -DBUILD_MODULE_Draw:BOOL=OFF \
      -DBUILD_MODULE_DETools:BOOL=OFF \
      -DUSE_FREETYPE:BOOL=OFF \
      -DUSE_OPENGL:BOOL=OFF \
      -DUSE_XLIB:BOOL=OFF \
      -DBUILD_DOC_Overview:BOOL=OFF \
      ${SUBPROJECT_CMAKE_ARGS}

emmake make -j9 install
