set -e
cd 
cd src/ngsolve
git submodule update --init --recursive
cd
mkdir -p build/ngsolve
cd build/ngsolve
cmake \
  ../../src/ngsolve \
  -DCMAKE_CXX_COMPILER=g++-7 \
  -DCMAKE_C_COMPILER=gcc-7 \
  -DCMAKE_CXX_FLAGS="-Og -Wall -Wno-sign-compare -DDebug" \
  -DUSE_NATIVE_ARCH=OFF \
  -DUSE_OCC=ON \
  -DUSE_CCACHE=ON \
  -DUSE_MKL=ON \
  -DUSE_UMFPACK=ON \
  -DINSTALL_PROFILES=ON \
  -DMKL_STATIC=ON \
  -DENABLE_UNIT_TESTS=ON \
  -DCMAKE_BUILD_TYPE=Debug \
  -DCPACK_PACKAGING_INSTALL_PREFIX=/opt/netgen \
  $CMAKE_ARGS
make -j12
make install
