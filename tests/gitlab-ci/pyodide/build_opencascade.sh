#! /bin/bash

set -e

source /root/emsdk/emsdk_env.sh

export EMSCRIPTEN_SYSROOT=$(em-config CACHE)/sysroot
export EMSCRIPTEN_INCLUDE=$EMSCRIPTEN_SYSROOT/include
export EMSCRIPTEN_BIN=$EMSCRIPTEN_SYSROOT/bin
export EMSCRIPTEN_LIB=$EMSCRIPTEN_SYSROOT/lib/wasm32-emscripten/pic
export TARGETINSTALLDIR=/root/xbuildenv/pyodide-root/cpython/installs/python-${PYTHON_VERSION}/
export CMAKE_TOOLCHAIN_FILE=/root/emsdk/upstream/emscripten/cmake/Modules/Platform/Emscripten.cmake

cd
wget https://github.com/Open-Cascade-SAS/OCCT/archive/refs/tags/V${OCC_VERSION}.zip
unzip V${OCC_VERSION}.zip
cd OCCT-${OCC_VERSION}
patch -p1 < /root/occ.patch
sed -i 's/-fexceptions/-fwasm-exceptions/g' adm/cmake/occt_defs_flags.cmake
mkdir build
cd build

emcmake cmake .. \
      -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_CROSSCOMPILING=ON \
      -DCMAKE_C_FLAGS="-I$TARGETINSTALLDIR/include/python3.13 -g0 -sSUPPORT_LONGJMP=0 -fwasm-exceptions -O2 -fvisibility=hidden" \
      -DCMAKE_CXX_FLAGS="-I$TARGETINSTALLDIR/include/python3.13 -g0 -sSUPPORT_LONGJMP=0 -fwasm-exceptions -O2 -fvisibility=hidden" \
      -DCMAKE_CXX_FLAGS_RELEASE="" \
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

set -o pipefail
emmake make -j9 install
# export VERBOSE=1
# emmake make -j9 install 2>&1 | tee /root/build_occ.log
cd /root/
rm -rf OCCT-${OCC_VERSION}
