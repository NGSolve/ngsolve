#!/bin/bash
export OSX_DEPLOYMENT_TARGET=10.15

PLATFORM=`python -c "import platform; print(platform.system())"`
if [ "$PLATFORM" = "Darwin" ]
then
	SHARED_EXT=dylib
else echo "not darwin"
	SHARED_EXT=so
fi

##############################################
# build netgen
mkdir build_netgen
cd build_netgen

cmake -G "Unix Makefiles" \
  -DUSE_CCACHE=ON \
  -DCMAKE_OSX_SYSROOT=/Library/Developer/CommandLineTools/SDKs/MacOSX10.15.sdk \
  -DCMAKE_OSX_DEPLOYMENT_TARGET=10.15 \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_PREFIX_PATH=${PREFIX} \
  -DCMAKE_INSTALL_PREFIX=${PREFIX} \
  -DCMAKE_SYSTEM_PREFIX_PATH=${PREFIX} \
  -DUSE_NATIVE_ARCH=OFF \
  -DUSE_SUPERBUILD=OFF \
  -DUSE_MPI=ON \
  -DUSE_GUI=ON \
  -DUSE_JPEG=ON \
  -DUSE_OCC=ON \
  -DUSE_MPEG=OFF \
  -DNG_INSTALL_DIR_BIN=bin \
  -DNG_INSTALL_DIR_LIB=lib/netgen \
  -DNG_INSTALL_DIR_CMAKE=lib/cmake/netgen \
  -DNG_INSTALL_DIR_PYTHON=${SP_DIR} \
  -DNG_INSTALL_DIR_RES=share \
  -DNG_INSTALL_DIR_INCLUDE=include/netgen \
  -DTCL_LIBRARY:FILEPATH=${PREFIX}/lib/libtcl8.6.${SHARED_EXT} \
  -DTK_LIBRARY:FILEPATH=${PREFIX}/lib/libtk8.6.${SHARED_EXT} \
  -DTCL_INCLUDE_PATH:FILEPATH=${SRC_DIR}/tcl/generic \
  -DTK_INCLUDE_PATH:FILEPATH=${SRC_DIR}/tk/generic \
  -DOCC_INCLUDE_DIR=${PREFIX}/include/opencascade \
  -DPYBIND11_INCLUDE_DIR:FILEPATH=${PREFIX}/include/python3.7m \
  -DPYTHON_INCLUDE_DIRS:FILEPATH=${PREFIX}/include/python3.7m \
  -DNG_INSTALL_PYBIND:BOOL=OFF \
  -DBUILD_STUB_FILES=OFF \
  -DBUILD_FOR_CONDA=ON \
  ${SRC_DIR}

make -j3
make install

