#!/bin/bash
export OSX_DEPLOYMENT_TARGET=10.12

PLATFORM=`python -c "import platform; print(platform.system())"`
if [ "$PLATFORM" = "Darwin" ]
then
	SHARED_EXT=dylib
else echo "not darwin"
	SHARED_EXT=so
fi

env
echo "BUILD_DIR" `pwd`
echo "SRC_DIR" `${SRC_DIR}`
export CCACHE_BASEDIR=$PREFIX/..

cd ..

##############################################
# build ngsolve
mkdir -p build_ngsolve
cd build_ngsolve

cmake -G "Unix Makefiles" \
  -DUSE_CCACHE=ON \
  -DUSE_SUPERBUILD=OFF \
  -DCMAKE_OSX_DEPLOYMENT_TARGET=10.12 \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX=${PREFIX} \
  -DUSE_NATIVE_ARCH=OFF \
  -DUSE_CCACHE=ON \
  -DUSE_MKL=ON \
  -DMKL_ROOT=${PREFIX} \
  -DMKL_SDL=ON \
  -DLAPACK_LIBRARIES=${PREFIX}/lib/libmkl_rt.${SHARED_EXT} \
  -DUSE_UMFPACK=OFF \
  -DBUILD_STUB_FILES=OFF \
   -DBUILD_JUPYTER_WIDGETS=ON \
  ${SRC_DIR}

make -j$CPU_COUNT
make install

